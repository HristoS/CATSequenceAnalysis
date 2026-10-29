using System.Globalization;
using System.Security.Claims;
using System.Text.Encodings.Web;
using Microsoft.AspNetCore.Authentication;
using Microsoft.Extensions.Options;
using Microsoft.OpenApi;
using Microsoft.AspNetCore.OpenApi;
using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Api
{
    /// <summary>
    /// D3: every request carries <c>Authorization: Bearer sp_...</c>, a key issued by the Ingest's <c>add-user</c>. Only the key's
    /// hash is stored (<see cref="Users"/>). The user id becomes the NameIdentifier claim, which every endpoint passes on as the
    /// catalogue's viewer.
    /// </summary>
    internal sealed class ApiKeyAuthentication(IOptionsMonitor<AuthenticationSchemeOptions> options, ILoggerFactory logger, UrlEncoder encoder,
        NpgsqlDataSource dataSource) : AuthenticationHandler<AuthenticationSchemeOptions>(options, logger, encoder)
    {
        public const string SchemeName = "ApiKey";

        private const string bearerPrefix = "Bearer ";

        protected override async Task<AuthenticateResult> HandleAuthenticateAsync()
        {
            var header = Request.Headers.Authorization.ToString();
            if (!header.StartsWith(bearerPrefix, StringComparison.OrdinalIgnoreCase)) return AuthenticateResult.NoResult();

            var user = await Users.AuthenticateAsync(dataSource, header[bearerPrefix.Length..].Trim(), Context.RequestAborted);
            if (user is null) return AuthenticateResult.Fail("Unknown or disabled API key.");

            var identity = new ClaimsIdentity(
                [new Claim(ClaimTypes.NameIdentifier, user.Id.ToString(CultureInfo.InvariantCulture)), new Claim(ClaimTypes.Name, user.Name)],
                SchemeName);
            return AuthenticateResult.Success(new AuthenticationTicket(new ClaimsPrincipal(identity), SchemeName));
        }

        protected override Task HandleChallengeAsync(AuthenticationProperties properties)
        {
            Response.StatusCode = StatusCodes.Status401Unauthorized;
            Response.Headers.WWWAuthenticate = "Bearer";
            return Task.CompletedTask;
        }
    }

    internal static class UserClaims
    {
        public static long UserId(this ClaimsPrincipal user) =>
            long.Parse(user.FindFirstValue(ClaimTypes.NameIdentifier)!, CultureInfo.InvariantCulture);
    }

    /// <summary>Declares the bearer key in the OpenAPI document, so clients generated from it send one.</summary>
    internal sealed class BearerSecurityTransformer : IOpenApiDocumentTransformer
    {
        public Task TransformAsync(OpenApiDocument document, OpenApiDocumentTransformerContext context, CancellationToken cancellationToken)
        {
            document.Components ??= new OpenApiComponents();
            document.Components.SecuritySchemes ??= new Dictionary<string, IOpenApiSecurityScheme>();
            document.Components.SecuritySchemes[ApiKeyAuthentication.SchemeName] = new OpenApiSecurityScheme
            {
                Type = SecuritySchemeType.Http,
                Scheme = "bearer",
                Description = "An API key from the Ingest's add-user command (sp_...).",
            };
            document.Security = [new OpenApiSecurityRequirement { [new OpenApiSecuritySchemeReference(ApiKeyAuthentication.SchemeName, document)] = [] }];
            return Task.CompletedTask;
        }
    }
}
