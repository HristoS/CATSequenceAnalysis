using Aspire.Hosting.Docker.Resources.ComposeNodes;
using Aspire.Hosting.Docker.Resources.ServiceNodes;
using Microsoft.Extensions.Configuration;
using Projects;

// The platform as a set of services (P6): `dotnet run --project SequencePlatform.AppHost` starts them with the dashboard, and the
// publish operation writes the same set as Docker Compose. Run mode uses the existing local database, through the connection string the
// Ingest and the tests already read; publish mode brings its own empty postgres:17.
var builder = DistributedApplication.CreateBuilder(args);

// No dashboard container: it would publish a port beyond this machine, and `docker compose logs` covers a local run.
builder.AddDockerComposeEnvironment("compose").WithDashboard(false);

builder.Configuration.AddJsonFile(Path.Combine(builder.AppHostDirectory, "..", "SequencePlatform.Ingest", "appsettings.json"), optional: true);

var server = builder.ExecutionContext.IsPublishMode
    ? builder.AddPostgres("postgres").WithImageTag("17").WithDataVolume()
        .PublishAsDockerComposeService((_, service) => service.Healthcheck = new Healthcheck
        {
            // Over TCP: the image's temporary server during first-time setup listens on a socket only, and migrate must not reach that one.
            Test = ["CMD-SHELL", "pg_isready -h 127.0.0.1 -U postgres"],
            Interval = "2s",
            Timeout = "3s",
            Retries = 30,
            StartPeriod = "2s",
        })
    : null;

IResourceBuilder<IResourceWithConnectionString> database = server is null
    ? builder.AddConnectionString("SequencePlatform")
    : server.AddDatabase("SequencePlatform", "sequence_platform");

// The API and the worker connect as sequence_app (migration 0021: reads everything, writes only uploads and search jobs), never as the
// owner, once hosted (S6.4). migrate, as the owner, gives that role a login with a generated password (persisted with the deployment's
// other secrets; letters and digits only, so it needs no quoting in a connection string). Run mode keeps the local connection string.
var appPassword = server is null
    ? null
    : builder.AddParameter("sequence-app-password", new GenerateParameterDefault { MinLength = 32, Special = false }, secret: true, persist: true);
// The role's name as Database.AppRole spells it (the AppHost doesn't reference the Core library).
const string appRole = "sequence_app";
IResourceBuilder<IResourceWithConnectionString> appDatabase = server is null
    ? database
    : builder.AddConnectionString("SequencePlatformApp", ReferenceExpression.Create(
        $"Host={server.Resource.PrimaryEndpoint.Property(EndpointProperty.Host)};Port={server.Resource.PrimaryEndpoint.Property(EndpointProperty.Port)};" +
        $"Username={appRole};Password={appPassword!.Resource};Database=sequence_platform"));

var migrate = builder.AddProject<SequencePlatform_Ingest>("migrate")
    .WithArgs("migrate")
    .WithEnvironment("SEQUENCE_PLATFORM_DB", database);
if (server is not null)
    migrate.WaitFor(server)
        .WithEnvironment("SEQUENCE_APP_PASSWORD", appPassword!)
        .PublishAsDockerComposeService((_, service) => service.DependsOn[server.Resource.Name] = new ServiceDependency { Condition = "service_healthy" });

// Localhost only: Caddy in front holds the TLS (hosting/tencent/Caddyfile).
builder.AddProject<SequencePlatform_Api>("api")
    .WithEnvironment("ConnectionStrings__SequencePlatform", appDatabase)
    .WaitForCompletion(migrate)
    .WithEnvironment("SearchWorker__Enabled", "false")
    .WithExternalHttpEndpoints()
    .PublishAsDockerComposeService((_, service) => service.Ports = ["127.0.0.1:5080:${API_PORT}"]);

// The same project with the worker on and no port anyone calls. Exactly one: at startup it requeues every running job as its own.
var worker = builder.AddProject<SequencePlatform_Api>("worker", launchProfileName: null)
    .WithEnvironment("ConnectionStrings__SequencePlatform", appDatabase)
    .WaitForCompletion(migrate)
    .WithEnvironment("ASPNETCORE_URLS", "http://127.0.0.1:0");
// Locally, localization keeps its levels in memory: 1.4-9.6x faster for about 0.8 GB (results/SearchBreakdown_Before.csv). The published
// set (the 8 GB server) keeps streaming.
if (builder.ExecutionContext.IsRunMode) worker.WithEnvironment("SearchWorker__CacheLocalization", "true");

builder.Build().Run();
