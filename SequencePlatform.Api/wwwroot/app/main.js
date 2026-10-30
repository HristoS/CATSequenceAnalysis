import { auth, get } from "./api.js";
import { CataloguePage } from "./catalogue.js";
import { Message, html, useEffect, useState } from "./html.js";
import { MapPage } from "./map.js";
import { RecordPage } from "./record.js";
import { SearchPage, SearchRunPage } from "./search.js";
import { TreePage } from "./tree.js";
import { UploadPage } from "./upload.js";

const links = [["catalogue", "Catalogue"], ["map", "Map"], ["tree", "Tree"], ["search", "Search"], ["upload", "Upload"]];

function parseRoute() {
  const [path, query = ""] = (location.hash.slice(2) || "catalogue").split("?");
  return { parts: path.split("/"), params: new URLSearchParams(query) };
}

function useRoute() {
  const [route, setRoute] = useState(parseRoute);
  useEffect(() => {
    const update = () => { setRoute(parseRoute()); window.scrollTo(0, 0); };
    window.addEventListener("hashchange", update);
    return () => window.removeEventListener("hashchange", update);
  }, []);
  return route;
}

function SignIn({ onSignedIn }) {
  const [key, setKey] = useState("");
  const [error, setError] = useState(null);

  async function submit(event) {
    event.preventDefault();
    setError(null);
    auth.key = key.trim();
    try {
      await get("/api/sequences?limit=1");
      onSignedIn();
    } catch (e) {
      auth.key = null;
      setError(e.message);
    }
  }

  return html`<main class="signin">
    <h1>Sequence Platform</h1>
    <p>Search and organize genomes from their Z-curve profiles. The platform stores profiles, never DNA.</p>
    <form onSubmit=${submit} class="stack">
      <label class="field"><span>API key</span>
        <input type="password" autoComplete="off" autoFocus value=${key} onInput=${e => setKey(e.target.value)} placeholder="sp_…" required />
      </label>
      <button class="primary" type="submit">Sign in</button>
    </form>
    ${error && html`<${Message} error>${error}<//>`}
    <p class="hint">Keys come from the Ingest's <code>add-user</code> command. This tab keeps yours until it closes.</p>
  </main>`;
}

function Page({ route }) {
  const [section, id] = route.parts;
  switch (section) {
    case "catalogue": return html`<${CataloguePage} />`;
    case "sequences": return html`<${RecordPage} key=${id} id=${id} />`;
    case "map": return html`<${MapPage} params=${route.params} />`;
    case "tree": return html`<${TreePage} params=${route.params} />`;
    case "search": return html`<${SearchPage} params=${route.params} />`;
    case "searches": return html`<${SearchRunPage} key=${id} id=${id} />`;
    case "upload": return html`<${UploadPage} />`;
    default: return html`<${Message} error>There is no page here.<//>`;
  }
}

function App() {
  const [signedIn, setSignedIn] = useState(!!auth.key);
  const route = useRoute();
  useEffect(() => auth.onSignedOut(() => setSignedIn(false)), []);

  if (!signedIn) return html`<${SignIn} onSignedIn=${() => setSignedIn(true)} />`;

  const section = route.parts[0];
  const active = section === "sequences" ? "catalogue" : section === "searches" ? "search" : section;
  return html`<div class="app">
    <header class="topbar">
      <a class="brand" href="#/catalogue">Sequence Platform</a>
      <nav aria-label="Main">
        ${links.map(([path, text]) => html`<a key=${path} href="#/${path}" aria-current=${active === path ? "page" : undefined}>${text}</a>`)}
      </nav>
      <button type="button" class="quiet" onClick=${() => { auth.key = null; setSignedIn(false); }}>Sign out</button>
    </header>
    <main class="page"><${Page} route=${route} /></main>
  </div>`;
}

window.ReactDOM.createRoot(document.getElementById("root")).render(html`<${App} />`);
