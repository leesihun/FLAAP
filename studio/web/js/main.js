// FLAAP Studio entry point: tab switching, hash routing, health badge.
//   #model                      -> model workspace
//   #results[/<run>[/<view>]]   -> results workspace (run = output/<run>)

import { api } from './api.js';
import { $, $$, toast } from './util.js';
import { initModelTab, modelTabShown } from './model_tab.js';
import { initResultsTab, loadRuns, openRun, showView, currentRun, currentView } from './results_tab.js';

let tab = 'model';
let runsLoaded = false;
let suppressHash = false;

function setHash(h) {
  if (location.hash === h) return;
  suppressHash = true;
  history.replaceState(null, '', h);
  suppressHash = false;
}

async function showTab(name, { run, view } = {}) {
  tab = name;
  for (const b of $$('#mainTabs button')) b.classList.toggle('active', b.dataset.tab === name);
  $('#ws-model').hidden = name !== 'model';
  $('#ws-results').hidden = name !== 'results';
  if (name === 'model') { modelTabShown(); setHash('#model'); return; }
  if (!runsLoaded || run) { runsLoaded = true; await loadRuns(); }
  const cur = currentRun();
  if (run && (!cur || cur.name !== run)) await openRun(run, view);
  else if (run && view) await showView(view);
  else if (!cur) {
    const first = $('#runList .item');
    if (first) first.click();
  }
  const c = currentRun();
  setHash(c ? `#results/${encodeURIComponent(c.name)}/${currentView()}` : '#results');
}

function parseHash() {
  const parts = location.hash.replace(/^#/, '').split('/').map(decodeURIComponent);
  return { tab: parts[0] === 'results' ? 'results' : 'model', run: parts[1] || undefined, view: parts[2] || undefined };
}

async function health() {
  const b = $('#serverBadge');
  try {
    const h = await api.health();
    b.textContent = `server ok · py ${h.python || ''}`.trim();
    b.className = 'badge ok';
    b.title = h.root || '';
  } catch {
    b.textContent = 'server offline';
    b.className = 'badge err';
  }
}

async function main() {
  for (const b of $$('#mainTabs button')) b.addEventListener('click', () => showTab(b.dataset.tab));
  initResultsTab({
    onRoute: (run, view) => { if (tab === 'results') setHash(`#results/${encodeURIComponent(run)}/${view}`); },
  });
  health();
  setInterval(health, 15000);
  const h = parseHash();
  const modelReady = initModelTab({
    openResults: async (name) => { runsLoaded = true; await loadRuns(); await showTab('results', { run: name, view: 'summary' }); },
    onRunFinished: () => { if (runsLoaded) loadRuns(); },
  }).catch((e) => toast(`model tab: ${e.message}`, true));
  if (h.tab === 'results') await showTab('results', h);
  else await modelReady;
  window.addEventListener('hashchange', () => {
    if (suppressHash) return;
    const q = parseHash();
    showTab(q.tab, q);
  });
}

main();
