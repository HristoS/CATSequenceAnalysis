using System.Net;
using Microsoft.AspNetCore.Mvc.Testing;

namespace SequencePlatform.Tests
{
    /// <summary>The web UI is static files in front of the API: open without a key, without weakening anything behind it.</summary>
    public class WebUiTests
    {
        private WebApplicationFactory<SequencePlatform.Api.SequencePageResponse> factory = null!;

        [OneTimeSetUp]
        public void Start() => factory = ApiHost.Create();

        [OneTimeTearDown]
        public void Stop() => factory.Dispose();

        [Test]
        public async Task ThePageLoadsWithoutAKeyAndLoadsNothingFromElsewhere()
        {
            var anonymous = factory.CreateClient();

            var page = await anonymous.GetAsync("/");
            Assert.That(page.StatusCode, Is.EqualTo(HttpStatusCode.OK));
            Assert.That(page.Content.Headers.ContentType!.MediaType, Is.EqualTo("text/html"));
            var body = await page.Content.ReadAsStringAsync();
            Assert.That(body, Does.Contain("default-src 'self'"), "the browser refuses anything from another origin");
            Assert.That(body, Does.Not.Contain("http://").And.Not.Contain("https://"));

            foreach (var path in new[] { "/app/main.js", "/app/experiments.js", "/app/zcurve3d.js", "/app/styles.css", "/vendor/react.production.min.js", "/vendor/react-dom.production.min.js", "/vendor/htm.mjs" })
                Assert.That((await anonymous.GetAsync(path)).StatusCode, Is.EqualTo(HttpStatusCode.OK), path);
        }

        [Test]
        public async Task EverythingBehindItStillNeedsAKey()
        {
            var anonymous = factory.CreateClient();

            foreach (var path in new[] { "/api/sequences", "/api/landscape", "/api/trees", "/api/no-such-route" })
                Assert.That((await anonymous.GetAsync(path)).StatusCode, Is.EqualTo(HttpStatusCode.Unauthorized), path);
        }

        [Test]
        public async Task HealthEndpointsAnswerWithoutAKeyAndSayNothingElse()
        {
            var anonymous = factory.CreateClient();

            foreach (var path in new[] { "/health", "/alive" })
                Assert.That(await anonymous.GetStringAsync(path), Is.EqualTo("Healthy"), path);
        }
    }
}
