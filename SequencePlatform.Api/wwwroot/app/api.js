// The key lives in sessionStorage: it goes with the tab, and it is never written anywhere else.
const keyName = "sequence-platform-key";
let signedOut = () => {};

export const auth = {
  get key() {
    try { return sessionStorage.getItem(keyName); } catch { return null; }
  },
  set key(value) {
    try {
      if (value) sessionStorage.setItem(keyName, value);
      else sessionStorage.removeItem(keyName);
    } catch { /* a blocked store only means signing in again */ }
  },
  onSignedOut(callback) { signedOut = callback; },
};

export const queryString = params => {
  const entries = Object.entries(params).filter(([, v]) => v !== "" && v != null);
  return entries.length ? `?${new URLSearchParams(entries)}` : "";
};

// Problem details carry a message written for the caller; anything else falls back to the status.
export async function send(path, { method = "GET", body, contentType } = {}) {
  const headers = { Authorization: `Bearer ${auth.key}` };
  if (contentType) headers["Content-Type"] = contentType;
  const response = await fetch(path, { method, body, headers });
  if (response.status === 401) {
    auth.key = null;
    signedOut();
    throw new Error("That key isn't valid.");
  }
  if (!response.ok) {
    const problem = await response.json().catch(() => null);
    throw new Error(problem?.detail ?? problem?.title ?? `The server answered ${response.status}.`);
  }
  return response.status === 204 ? null : response.json();
}

export const get = path => send(path);
