// Thin wrapper around the studio JSON API.

async function req(url, opts) {
  const r = await fetch(url, opts);
  const ct = r.headers.get('content-type') || '';
  const body = ct.includes('json') ? await r.json() : await r.text();
  if (!r.ok) throw new Error((body && body.error) || `${r.status} ${r.statusText}`);
  return body;
}
const post = (url, data) => req(url, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(data) });

export const fileUrl = (path) => `/api/file?path=${encodeURIComponent(path)}`;

const cache = new Map();
export const api = {
  health: () => req('/api/health'),
  models: () => req('/api/models'),
  presets: () => req('/api/presets'),
  configs: () => req('/api/configs'),
  config: (path) => req(`/api/config?path=${encodeURIComponent(path)}`),
  validate: (text, source) => post('/api/validate', { text, source }),
  run: (text, name, source) => post('/api/run', { text, name, source }),
  job: (id, since = 0) => req(`/api/job?id=${encodeURIComponent(id)}&since=${since}`),
  cancel: (id) => post('/api/cancel', { id }),
  outputs: () => req('/api/outputs'),
  /** JSON from output/ (cached per path + mtime tag). */
  json: async (path, tag = '') => {
    const key = path + '|' + tag;
    if (!cache.has(key)) cache.set(key, req(fileUrl(path)).catch((e) => { cache.delete(key); throw e; }));
    return cache.get(key);
  },
  /** A FLAAP result text file (results.txt, modes.txt, flutter_*.txt ...) parsed by the server. */
  result: async (path, tag = '') => {
    const key = 'R|' + path + '|' + tag;
    if (!cache.has(key)) cache.set(key, req(`/api/result?path=${encodeURIComponent(path)}`).catch((e) => { cache.delete(key); throw e; }));
    return cache.get(key);
  },
  text: (path) => req(fileUrl(path)),
};
