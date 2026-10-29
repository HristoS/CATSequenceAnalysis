using Microsoft.AspNetCore.Diagnostics;
using Microsoft.AspNetCore.Mvc;

namespace SequencePlatform.Api
{
    /// <summary>
    /// Core validates its inputs with ArgumentException (an unknown rank, a page size out of range, a profile whose declared length
    /// disagrees with its windows) and rejects malformed FASTA or ZCWP with FormatException. Those messages are written for the caller
    /// and carry no sequence, so they become a 400 with the message as detail. Kestrel's own request errors (a body over the size
    /// limit: 413) keep their status. Anything else is a 500 from the default handler, without detail.
    /// </summary>
    internal sealed class ApiExceptionHandler(IProblemDetailsService problemDetails, ILogger<ApiExceptionHandler> logger) : IExceptionHandler
    {
        public async ValueTask<bool> TryHandleAsync(HttpContext httpContext, Exception exception, CancellationToken cancellationToken)
        {
            var (status, title, detail) = exception switch
            {
                // Without the " (Parameter 'filter')" suffix, which names a C# parameter, not a query parameter.
                ArgumentException a => (StatusCodes.Status400BadRequest, "Invalid request",
                    a.ParamName is null ? a.Message : a.Message.Replace($" (Parameter '{a.ParamName}')", "")),
                FormatException f => (StatusCodes.Status400BadRequest, "Invalid upload", f.Message),
                BadHttpRequestException b => (b.StatusCode, "Bad request", b.Message),
                _ => (0, null, null),
            };
            if (status == 0) return false;

            logger.LogInformation("Rejected request {Path} with {Status}: {Detail}", httpContext.Request.Path, status, detail);
            httpContext.Response.StatusCode = status;
            return await problemDetails.TryWriteAsync(new ProblemDetailsContext
            {
                HttpContext = httpContext,
                Exception = exception,
                ProblemDetails = new ProblemDetails { Status = status, Title = title, Detail = detail },
            });
        }
    }
}
