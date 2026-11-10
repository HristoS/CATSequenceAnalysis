using System.Text.Json.Serialization;
using Microsoft.AspNetCore.Authentication;
using Microsoft.AspNetCore.Authorization;
using Npgsql;
using SequencePlatform.Api;
using SequencePlatform.Core;

// The platform's HTTP API (P4). Kestrel serves plain HTTP; hosted, it listens on localhost behind Caddy's TLS and connects to PostgreSQL
// as the least-privilege role sequence_app (S6.4, migration 0021). The connection string comes from ASP.NET configuration (user-secrets, or the
// ConnectionStrings__SequencePlatform variable), else SEQUENCE_PLATFORM_DB as the Ingest reads it. The API doesn't migrate;
// run the Ingest's migrate first.
var builder = WebApplication.CreateBuilder(args);
builder.AddServiceDefaults();

// Uploads stream, but a request still has an upper bound: 32 MiB holds the largest known prokaryote (~13.7 Mbp) as FASTA.
builder.WebHost.ConfigureKestrel(kestrel => kestrel.Limits.MaxRequestBodySize = 32 * 1024 * 1024);

builder.Services.AddSingleton(services =>
{
    var connectionString = services.GetRequiredService<IConfiguration>().GetConnectionString("SequencePlatform") ?? ConnectionString.Resolve();
    return NpgsqlDataSource.Create(Database.ConnectionSettings(connectionString).ConnectionString);
});

builder.Services.AddAuthentication(ApiKeyAuthentication.SchemeName)
    .AddScheme<AuthenticationSchemeOptions, ApiKeyAuthentication>(ApiKeyAuthentication.SchemeName, null);
builder.Services.AddAuthorization(options =>
    options.FallbackPolicy = new AuthorizationPolicyBuilder().RequireAuthenticatedUser().Build());

builder.Services.AddSingleton<SignatureIndex>();
builder.Services.AddSingleton<Landscape>();
builder.Services.Configure<SearchWorkerOptions>(builder.Configuration.GetSection("SearchWorker"));
builder.Services.AddHostedService<SearchWorker>();

// The map of everything is several MB of JSON. Plain HTTP only: behind a TLS proxy the request counts as HTTPS (forwarded headers), so the
// proxy compresses instead (hosting/tencent/Caddyfile, encode).
builder.Services.AddResponseCompression();
builder.Services.AddProblemDetails();
builder.Services.AddExceptionHandler<ApiExceptionHandler>();
builder.Services.AddOpenApi(options => options.AddDocumentTransformer<BearerSecurityTransformer>());
builder.Services.ConfigureHttpJsonOptions(options =>
{
    options.SerializerOptions.Converters.Add(new JsonStringEnumConverter());
    options.SerializerOptions.DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull;
    options.SerializerOptions.NumberHandling = JsonNumberHandling.Strict;
    options.SerializerOptions.AllowDuplicateProperties = false;
});

var app = builder.Build();
_ = app.Services.GetRequiredService<SignatureIndex>().RefreshAsync();   // warm the relatives index without delaying startup

app.UseExceptionHandler();
app.UseStatusCodePages();
app.UseResponseCompression();
app.UseDefaultFiles();
app.UseStaticFiles();   // the web UI, before authentication: it holds no data, and every call it makes carries the user's key
app.UseAuthentication();
app.UseAuthorization();

if (app.Environment.IsDevelopment()) app.MapOpenApi().AllowAnonymous();
app.MapDefaultEndpoints();
app.MapSequences();
app.MapSearches();
app.MapOrganisms();
app.MapGenes();
app.MapExperiments();

app.Run();
