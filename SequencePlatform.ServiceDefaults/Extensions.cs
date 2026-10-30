using Microsoft.AspNetCore.Builder;
using Microsoft.AspNetCore.Diagnostics.HealthChecks;
using Microsoft.AspNetCore.Routing;
using Microsoft.Extensions.DependencyInjection;
using Microsoft.Extensions.Diagnostics.HealthChecks;
using Microsoft.Extensions.Logging;
using OpenTelemetry;
using OpenTelemetry.Metrics;
using OpenTelemetry.Trace;

namespace Microsoft.Extensions.Hosting
{
    /// <summary>
    /// What every platform service shares under Aspire: logs, metrics and traces to the dashboard, and health endpoints. Telemetry leaves
    /// the process only when OTEL_EXPORTER_OTLP_ENDPOINT is set, which the AppHost does. Database calls are deliberately not traced, and
    /// ASP.NET's instrumentation drops query strings, so a span can't carry a record's name or anything else a request sent.
    /// </summary>
    public static class Extensions
    {
        private const string healthPath = "/health";
        private const string alivePath = "/alive";

        public static TBuilder AddServiceDefaults<TBuilder>(this TBuilder builder) where TBuilder : IHostApplicationBuilder
        {
            builder.Logging.AddOpenTelemetry(logging =>
            {
                logging.IncludeFormattedMessage = true;
                logging.IncludeScopes = true;
            });

            builder.Services.AddOpenTelemetry()
                .WithMetrics(metrics => metrics.AddAspNetCoreInstrumentation().AddHttpClientInstrumentation().AddRuntimeInstrumentation())
                .WithTracing(tracing => tracing
                    .AddSource(builder.Environment.ApplicationName)
                    .AddAspNetCoreInstrumentation(options => options.Filter = context =>
                        !context.Request.Path.StartsWithSegments(healthPath) && !context.Request.Path.StartsWithSegments(alivePath))
                    .AddHttpClientInstrumentation());

            if (!string.IsNullOrWhiteSpace(builder.Configuration["OTEL_EXPORTER_OTLP_ENDPOINT"]))
                builder.Services.AddOpenTelemetry().UseOtlpExporter();

            builder.Services.AddHealthChecks().AddCheck("self", () => HealthCheckResult.Healthy(), ["live"]);
            return builder;
        }

        /// <summary>Open without a key: an orchestrator has none, and they say only whether the process is up.</summary>
        public static IEndpointRouteBuilder MapDefaultEndpoints(this IEndpointRouteBuilder app)
        {
            app.MapHealthChecks(healthPath).AllowAnonymous();
            app.MapHealthChecks(alivePath, new HealthCheckOptions { Predicate = check => check.Tags.Contains("live") }).AllowAnonymous();
            return app;
        }
    }
}
