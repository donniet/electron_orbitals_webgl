// Recording UI: the red dot over the view, and the recording controls that
// replace the periodic-table strip while they are open.
//
// States: closed -> setup (dot blinks, outside the crop dimmed)
//         -> rendering (dot hidden, live view frozen, progress bar)
//         -> done (download link) or error; back to setup with "Record again".

import {
  RESOLUTIONS, FPS_CHOICES, QUALITY, evenSize, paneCrops, rotationPeriod, framePlan,
  supportedFormats, estimateBytes, record, videoSupported,
} from './record.js';

const $ = sel => document.querySelector(sel);

function h(tag, attrs = {}, ...children) {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k.startsWith('on')) e.addEventListener(k.slice(2), v);
    else if (k === 'className') e.className = v;
    else if (v !== false && v != null) e.setAttribute(k, v === true ? '' : v);
  }
  for (const c of children.flat()) if (c != null) e.append(c);
  return e;
}

function formatBytes(b) {
  if (b >= 1e9) return `${(b / 1e9).toFixed(1)} GB`;
  if (b >= 1e6) return `${(b / 1e6).toFixed(b >= 1e8 ? 0 : 1)} MB`;
  return `${Math.max(1, Math.round(b / 1e3))} kB`;
}

// Above this estimate, warn that the browser may run out of memory.
const LARGE_BYTES = 1.5e9;

/**
 * ctx: {
 *   app,                          the app (state, tables)
 *   viewInfo() -> { width, height, bohrPerPixel, dpr, panes }   the on-screen view
 *                                 now; panes: pane rects (GL origin) as drawn
 *   views() -> [{ element, visible }]                    panes to record
 *   changed()                     redraw (crop preview, frozen view)
 * }
 * Returns { phase(), crop(), drawCrop(overlayCanvas), refresh(), open() }.
 */
export function createRecorderUI(ctx) {
  const { app } = ctx;
  const dot = $('#rec-dot');
  const bar = $('#recbar');
  const ptable = $('#ptable');

  const settings = { res: 'current', durMode: 'rotation', seconds: 5, fps: 60, grain: 2, quality: 'high', format: null };
  let phase = 'closed';                 // closed | setup | rendering | done | error
  let formats = [];
  let abort = null, result = null, url = null, error = null, progress = [0, 0];

  if (!videoSupported()) {
    dot.disabled = true;
    dot.title = 'Video recording needs WebCodecs (VideoEncoder), which this browser does not have.';
  }
  dot.addEventListener('click', () => (phase === 'closed' ? open() : phase !== 'rendering' && close()));

  function outSize() {
    if (settings.res === 'current') {
      const v = ctx.viewInfo();
      return evenSize(v.width, v.height);
    }
    const r = RESOLUTIONS.find(x => x.id === settings.res);
    return [r.w, r.h];
  }

  function duration() {
    return settings.durMode === 'rotation' && rotationPeriod(app.state.spin)
      ? { rotation: true } : { seconds: settings.seconds };
  }

  function frameCount() {
    return framePlan(app.state.orientation, app.state.spin, settings.fps, duration()).n;
  }

  function open() { phase = 'setup'; refreshFormats(); render(); }
  function close() {
    if (url) { URL.revokeObjectURL(url); url = null; }
    result = null; error = null; phase = 'closed'; render();
  }

  // Probe the encoders for the current size / fps / quality; keep the chosen
  // format if it is still available.
  let probe = 0;
  async function refreshFormats() {
    const id = ++probe;
    const [w, h] = outSize();
    const list = await supportedFormats(w, h, settings.fps, settings.quality);
    if (id !== probe) return;
    formats = list;
    if (!formats.some(f => f.id === settings.format)) settings.format = formats[0]?.id ?? null;
    if (phase === 'setup') render();
  }

  async function start() {
    const [w, h] = outSize();
    abort = new AbortController();
    phase = 'rendering';
    progress = [0, frameCount()];
    render();
    try {
      result = await record({
        tables: app.tables, state: { ...app.state }, views: ctx.views(), view: ctx.viewInfo(),
        width: w, height: h, fps: settings.fps, duration: duration(), grain: settings.grain,
        formatId: settings.format, quality: settings.quality, signal: abort.signal,
        onProgress: (i, n) => { progress = [i, n]; updateProgress(); },
      });
      url = URL.createObjectURL(result.blob);
      phase = 'done';
    } catch (e) {
      if (e.name === 'AbortError') phase = 'setup';
      else { error = e.message || String(e); phase = 'error'; }
    }
    abort = null;
    render();
  }

  // ---------------------------------------------------------------- DOM
  function field(label, control) {
    return h('label', { className: 'rb-field' }, h('span', {}, label), control);
  }

  function select(options, value, onchange) {
    const s = h('select', { onchange: e => onchange(e.target.value) },
      options.map(([v, text]) => h('option', { value: v, selected: String(v) === String(value) }, text)));
    return s;
  }

  function setupControls() {
    const v = ctx.viewInfo();
    const [cw, ch] = evenSize(v.width, v.height);
    const period = rotationPeriod(app.state.spin);
    const [w, hgt] = outSize();
    const n = frameCount();
    const bytes = estimateBytes(w, hgt, settings.fps, n, settings.quality);
    const set = (k, parse = x => x) => val => { settings[k] = parse(val); refreshFormats(); render(); };

    const durRotation = h('label', { title: period ? 'Exactly one turn: the video loops seamlessly.' : 'The atom is not spinning (flick it, or Reset view).' },
      h('input', { type: 'radio', name: 'rb-dur', checked: settings.durMode === 'rotation' && !!period, disabled: !period,
        onchange: () => { settings.durMode = 'rotation'; render(); } }),
      period ? ` 1 rotation (${period.toFixed(1)} s)` : ' 1 rotation');
    const secs = h('input', { type: 'number', min: 0.5, max: 600, step: 0.5, value: settings.seconds, className: 'rb-num',
      onchange: e => { settings.seconds = Math.min(600, Math.max(0.5, +e.target.value || 5)); settings.durMode = 'seconds'; render(); } });
    const durSeconds = h('label', {},
      h('input', { type: 'radio', name: 'rb-dur', checked: settings.durMode === 'seconds' || !period,
        onchange: () => { settings.durMode = 'seconds'; render(); } }), ' ', secs, ' s');

    const grain = h('span', { className: 'rb-radios', title: 'Output pixels per noise sample. 2 px keeps each detection\'s colour exact in compressed video and makes smaller files.' },
      [1, 2].map(g => h('label', {}, h('input', { type: 'radio', name: 'rb-grain', checked: settings.grain === g,
        onchange: () => { settings.grain = g; render(); } }), ` ${g} px`)));

    const formatSel = formats.length
      ? select(formats.map(f => [f.id, f.label]), settings.format, val => { settings.format = val; render(); })
      : h('span', { className: 'rb-muted' }, 'none available at this size');

    return [
      h('span', { className: 'rb-title' }, 'Record video'),
      field('Resolution', select([['current', `Current view (${cw} × ${ch})`], ...RESOLUTIONS.map(r => [r.id, r.label])],
        settings.res, set('res'))),
      field('Duration', h('span', { className: 'rb-radios' }, durRotation, durSeconds)),
      field('FPS', select(FPS_CHOICES.map(f => [f, String(f)]), settings.fps, set('fps', Number))),
      field('Grain', grain),
      field('Quality', select(Object.entries(QUALITY).map(([k, q]) => [k, q.label]), settings.quality, set('quality'))),
      field('Format', formatSel),
      h('span', { className: 'rb-muted' + (bytes > LARGE_BYTES ? ' rb-warn' : '') },
        `${n} frames · ≈ ${formatBytes(bytes)}` + (bytes > LARGE_BYTES ? ' (may not fit in memory)' : '')),
      h('span', { className: 'rb-actions' },
        h('button', { className: 'rb-primary', disabled: !settings.format, onclick: start, title: 'Render every frame and encode the video' }, '▶ Record'),
        h('button', { onclick: close, title: 'Close the recording controls' }, '✕')),
    ];
  }

  function updateProgress() {
    const p = bar.querySelector('progress'), t = bar.querySelector('.rb-progress-text');
    if (p) { p.max = progress[1]; p.value = progress[0]; }
    if (t) t.textContent = `Rendering frame ${progress[0]} / ${progress[1]}`;
  }

  function render() {
    dot.hidden = phase === 'rendering';
    dot.classList.toggle('blink', phase === 'setup' || phase === 'done' || phase === 'error');
    dot.setAttribute('aria-pressed', String(phase !== 'closed'));
    if (videoSupported()) dot.title = phase === 'closed' ? 'Record a video' : 'Close the recording controls';
    ptable.hidden = phase !== 'closed';
    bar.hidden = phase === 'closed';
    bar.replaceChildren();
    if (phase === 'setup') bar.append(...setupControls());
    else if (phase === 'rendering') {
      bar.append(h('span', { className: 'rb-title' }, 'Recording'),
        h('progress', { max: progress[1] || 1, value: progress[0] }),
        h('span', { className: 'rb-progress-text rb-muted' }),
        h('span', { className: 'rb-actions' }, h('button', { onclick: () => abort?.abort() }, 'Cancel')));
      updateProgress();
    } else if (phase === 'done') {
      bar.append(h('span', { className: 'rb-title' }, 'Video ready'),
        h('a', { className: 'rb-primary rb-download', href: url, download: result.filename },
          `⬇ Download ${result.filename} (${formatBytes(result.bytes)})`),
        h('span', { className: 'rb-actions' },
          h('button', { onclick: () => { URL.revokeObjectURL(url); url = null; result = null; open(); } }, 'Record again'),
          h('button', { onclick: close }, 'Close')));
    } else if (phase === 'error') {
      bar.append(h('span', { className: 'rb-title' }, 'Recording failed'),
        h('span', { className: 'rb-warn' }, error),
        h('span', { className: 'rb-actions' }, h('button', { onclick: open }, 'Back'), h('button', { onclick: close }, 'Close')));
    }
    ctx.changed();
  }

  // What the video will show: one rectangle per pane (device px, top-left
  // origin) while setting up, else null.
  function crop() {
    if (phase !== 'setup') return null;
    const v = ctx.viewInfo();
    const [w, hgt] = outSize();
    return paneCrops(v.width, v.height, v.panes, w, hgt, settings.grain);
  }

  // Dim everything outside the crop rectangles (live overlay only; never
  // part of the video). No outline: the blinking dot marks recording mode,
  // and with "Current view" nothing is dimmed at all.
  function drawCrop(canvas) {
    const rects = crop();
    if (!rects) return;
    const g = canvas.getContext('2d');
    g.save();
    g.fillStyle = app.state.mode === 'white' ? 'rgba(255,255,255,0.6)' : 'rgba(0,0,0,0.55)';
    g.beginPath();
    g.rect(0, 0, canvas.width, canvas.height);
    for (const c of rects) g.rect(c.x, c.y, c.w, c.h);
    g.fill('evenodd');
    g.restore();
  }

  // The view size or spin can change while the controls are open (resize,
  // device rotation, a fling): refresh sizes, formats and frame counts.
  function refresh() { if (phase === 'setup') { refreshFormats(); render(); } }

  return { phase: () => phase, crop, drawCrop, refresh, open: () => { if (phase === 'closed' && videoSupported()) open(); } };
}
