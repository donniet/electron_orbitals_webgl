// Video tests: encoders, muxers and recording. Media playback needs real
// time, so this page must NOT run under --virtual-time-budget. Serve the repo
// root and open http://localhost:8000/tests/video.html, or headless:
//   google-chrome --headless=new --no-sandbox --use-angle=swiftshader \
//     --enable-unsafe-swiftshader --user-data-dir="$(mktemp -d)" \
//     --enable-logging=stderr --v=0 http://localhost:8000/tests/video.html 2>&1 | grep -m1 -o 'TESTS.*'
// document.title / the console line become "PASS n/n" or "FAIL k/n".

import { muxMp4 } from '../js/mux-mp4.js';
import { muxWebm } from '../js/mux-webm.js';
import { cropRect, paneCrops, framePlan, rotationPeriod, h264Codec, vp9Codec, av1Codec, evenSize, record,
  supportedFormats, FORMATS } from '../js/record.js';
import { loadTables } from '../js/tables.js';
import { ELEMENTS, ORBITALS } from '../js/elements.js';
import { viewGeometry } from '../js/renderer.js';
import { qAxisAngle } from '../js/math.js';

const TABLES = await loadTables();

const results = [];
async function atest(name, fn) {
  try {
    await fn();
    results.push({ name, ok: true });
  } catch (e) {
    results.push({ name, ok: false, msg: e.message });
  }
}
function assert(cond, msg = 'assertion failed') { if (!cond) throw new Error(msg); }
function near(a, b, tol, msg = '') {
  if (!(Math.abs(a - b) <= tol)) throw new Error(`${msg} expected ${b} got ${a} (tol ${tol})`);
}

// ------------------------------------------------------------ video muxers
// Encode `n` solid-color frames with WebCodecs, mux them, then play the file
// in a <video>: size, duration and the decoded color at several frames.
const MUX_COLORS = [[230, 40, 40], [40, 200, 60], [40, 70, 230], [230, 220, 40], [220, 50, 220], [40, 220, 220]];

async function encodeColorFrames(codec, { width = 192, height = 108, fps = 30, n = 12, keyEvery = 4 } = {}) {
  const canvas = new OffscreenCanvas(width, height);
  const ctx = canvas.getContext('2d');
  const chunks = [];
  let description = null;
  const encoder = new VideoEncoder({
    output: (chunk, meta) => {
      const data = new Uint8Array(chunk.byteLength);
      chunk.copyTo(data);
      chunks.push({ data, pts: chunk.timestamp, key: chunk.type === 'key' });
      if (meta?.decoderConfig?.description) description = meta.decoderConfig.description;
    },
    error: e => { throw e; },
  });
  const config = { codec, width, height, framerate: fps, bitrate: 2e6 };
  if (codec.startsWith('avc1')) config.avc = { format: 'avc' };
  encoder.configure(config);
  for (let i = 0; i < n; i++) {
    const [r, g, b] = MUX_COLORS[i % MUX_COLORS.length];
    ctx.fillStyle = `rgb(${r},${g},${b})`;
    ctx.fillRect(0, 0, width, height);
    const frame = new VideoFrame(canvas, { timestamp: Math.round((i * 1e6) / fps), duration: Math.round(1e6 / fps) });
    encoder.encode(frame, { keyFrame: i % keyEvery === 0 });
    frame.close();
  }
  await encoder.flush();
  encoder.close();
  return { chunks, description, width, height, fps, n };
}

function videoEvent(video, type) {
  return new Promise((resolve, reject) => {
    const t = setTimeout(() => reject(new Error(`video: no ${type} event`)), 20000);
    video.addEventListener(type, () => { clearTimeout(t); resolve(); }, { once: true });
    video.addEventListener('error', () => { clearTimeout(t); reject(new Error(`video error ${video.error?.code} ${video.error?.message}`)); }, { once: true });
  });
}

async function checkPlayback(bytes, mime, { width, height, fps, n }) {
  const url = URL.createObjectURL(new Blob([bytes], { type: mime }));
  const video = document.createElement('video');
  video.muted = true;
  video.preload = 'auto';
  video.src = url;
  try {
    await videoEvent(video, 'loadeddata');
    assert(video.videoWidth === width && video.videoHeight === height, `size ${video.videoWidth}x${video.videoHeight}`);
    near(video.duration, n / fps, 1.5 / fps, 'duration');
    const c = document.createElement('canvas');
    c.width = width; c.height = height;
    const ctx = c.getContext('2d', { willReadFrequently: true });
    for (const i of [0, 3, 5, 8, n - 1]) {
      video.currentTime = (i + 0.5) / fps;
      await videoEvent(video, 'seeked');
      ctx.drawImage(video, 0, 0);
      const px = ctx.getImageData(width >> 1, height >> 1, 1, 1).data;
      const want = MUX_COLORS[i % MUX_COLORS.length];
      const err = Math.max(...want.map((v, j) => Math.abs(v - px[j])));
      assert(err < 24, `frame ${i}: got ${[...px.slice(0, 3)]} want ${want}`);
    }
  } finally {
    URL.revokeObjectURL(url);
  }
}

for (const [label, codec] of [['H.264', 'avc1.42E01F']]) {
  await atest(`MP4 (${label}): muxed file plays with the right size, duration and frames`, async () => {
    const e = await encodeColorFrames(codec);
    assert(e.description, 'encoder gave no avcC description');
    const timescale = e.fps * 1000;
    const bytes = muxMp4({ width: e.width, height: e.height, timescale, frameDuration: 1000, description: e.description,
      chunks: e.chunks.map(c => ({ ...c, pts: Math.round((c.pts * timescale) / 1e6) })) });
    await checkPlayback(bytes, 'video/mp4', e);
  });
}

for (const [label, codec] of [['VP9', 'vp09.00.10.08'], ['VP8', 'vp8'], ['AV1', 'av01.0.04M.08']]) {
  await atest(`WebM (${label}): muxed file plays with the right size, duration and frames`, async () => {
    const e = await encodeColorFrames(codec);
    const bytes = muxWebm({ width: e.width, height: e.height, fps: e.fps, codec, description: e.description, chunks: e.chunks });
    await checkPlayback(bytes, 'video/webm', e);
  });
}

// ------------------------------------------------------------ recording
await atest('cropRect: largest centred rectangle of the output aspect', async () => {
  const c = cropRect(1600, 900, 1080, 1080);
  near(c.w, 900, 1e-9); near(c.h, 900, 1e-9); near(c.x, 350, 1e-9); near(c.y, 0, 1e-9); near(c.scale, 900 / 1080, 1e-12);
  const d = cropRect(800, 1000, 1920, 1080);
  near(d.w, 800, 1e-9); near(d.h, 450, 1e-9); near(d.y, 275, 1e-9);
  assert(evenSize(1367, 769).join() === '1366,768', 'even size');
});

await atest('paneCrops: one centred crop per output pane, at the recording scale', async () => {
  const one = paneCrops(1600, 900, [{ x: 0, y: 0, w: 1600, h: 900 }], 1080, 1080);
  const c = cropRect(1600, 900, 1080, 1080);
  near(one[0].x, c.x, 1e-9); near(one[0].w, c.w, 1e-9); near(one[0].h, c.h, 1e-9);
  // Two panes side by side on screen, portrait output: panes stack in the
  // video, each crop centred on its own on-screen pane.
  const panes = [{ x: 0, y: 0, w: 800, h: 900 }, { x: 800, y: 0, w: 800, h: 900 }];
  const two = paneCrops(1600, 900, panes, 1080, 1920);
  const s = 900 / 1920;
  near(two[0].x + two[0].w / 2, 400, 1e-9); near(two[1].x + two[1].w / 2, 1200, 1e-9);
  near(two[0].y + two[0].h / 2, 450, 1e-9);
  near(two[0].w, 1080 * s, 1e-9); near(two[0].h, 960 * s, 1e-9);
});

await atest('framePlan: one rotation loops exactly; seconds use the spin rate', async () => {
  const q0 = [0.1, 0.2, 0.3, 0.9];
  const spin = { axis: [0, 1, 0], rate: 0.0025 };
  near(rotationPeriod(spin), 2 * Math.PI / 2.5, 1e-12);
  const p = framePlan(q0, spin, 60, { rotation: true });
  assert(p.n === Math.round(rotationPeriod(spin) * 60), `n ${p.n}`);
  const a = p.orientation(0), b = p.orientation(p.n);
  const same = Math.abs(a.reduce((s, v, i) => s + v * b[i], 0));   // q and -q are the same rotation
  near(same, 1, 1e-9, 'frame n == frame 0');
  const s = framePlan(q0, spin, 30, { seconds: 2 });
  assert(s.n === 60, `n ${s.n}`);
  near(s.step, 0.0025 * 1000 / 30, 1e-12);
  const still = framePlan(q0, null, 30, { seconds: 1 });
  assert(still.step === 0 && still.n === 30, 'no spin: still');
  let threw = false;
  try { framePlan(q0, null, 30, { rotation: true }); } catch { threw = true; }
  assert(threw, 'one rotation without spin must fail');
});

await atest('codec strings: levels by size, frame rate and bitrate', async () => {
  const eq = (got, want) => assert(got === want, `${got} != ${want}`);
  eq(h264Codec(1280, 720, 30), 'avc1.64001F');
  eq(h264Codec(1920, 1080, 30), 'avc1.640028');
  eq(h264Codec(1920, 1080, 60), 'avc1.64002A');
  eq(h264Codec(3840, 2160, 30), 'avc1.640033');
  eq(h264Codec(3840, 2160, 60), 'avc1.640034');
  eq(h264Codec(1920, 1080, 60, 100e6), 'avc1.640032');    // 100 Mbit/s needs level 5.0
  eq(vp9Codec(1920, 1080, 60), 'vp09.00.41.08');
  eq(vp9Codec(3840, 2160, 30), 'vp09.00.50.08');
  eq(vp9Codec(1920, 1080, 60, 100e6), 'vp09.00.51.08');
  eq(av1Codec(1920, 1080, 30), 'av01.0.08M.08');
  eq(av1Codec(3840, 2160, 60), 'av01.0.13M.08');
  eq(av1Codec(1920, 1080, 60, 45e6), 'av01.0.09H.08');     // over main tier 4.1, within high
});

await atest('supportedFormats: H.264 MP4 and VP9 WebM available at 1080p60', async () => {
  const ids = (await supportedFormats(1920, 1080, 60, 'high')).map(f => f.id);
  assert(ids.includes('mp4-h264') && ids.includes('webm-vp9'), ids.join());
});

// A short, deterministic recording of hydrogen spinning fast (12 frames per
// turn at 24 fps): output plays with the planned size, length and frames,
// and the source frames are centred (the crop keeps the atom in the middle).
function recState(over = {}) {
  return {
    element: ELEMENTS[0], orientation: [0, 0, 0, 1], zoom: 1, density: 1, equalize: false, kMax: 4,
    gamma: 0.5, mode: 'dark', fixedScale: false, visible: ORBITALS.map(() => true), cut: false,
    color: 'blend', spin: { axis: [0, 1, 0], rate: (2 * Math.PI) / 500 }, overlay: false, debugSeed: true, ...over,
  };
}
const VIEW = (() => {
  const g = viewGeometry(ELEMENTS[0], 400, 300, 1);
  return { width: 400, height: 300, bohrPerPixel: g.bohrPerPixel, dpr: 1 };
})();

for (const formatId of ['mp4-h264', 'webm-vp9']) {
  await atest(`record (${formatId}): one rotation, plays with the planned frames`, async () => {
    const state = recState();
    let progress = 0, centroids = [];
    const res = await record({
      tables: TABLES, state, views: [{ element: state.element, visible: state.visible }], view: VIEW,
      width: 320, height: 240, fps: 24, duration: { rotation: true }, grain: 2, formatId, quality: 'high',
      onProgress: i => { progress = i; },
      onFrame: (canvas, i) => {
        if (i > 1) return;
        const c = document.createElement('canvas');
        c.width = canvas.width; c.height = canvas.height;
        const ctx = c.getContext('2d');
        ctx.drawImage(canvas, 0, 0);
        const px = ctx.getImageData(0, 0, c.width, c.height).data;
        let sx = 0, sy = 0, sw = 0;
        for (let y = 0; y < c.height; y++) for (let x = 0; x < c.width; x++) {
          const v = px[4 * (y * c.width + x)] + px[4 * (y * c.width + x) + 1] + px[4 * (y * c.width + x) + 2];
          sx += v * (x + 0.5); sy += v * (y + 0.5); sw += v;
        }
        centroids.push([sx / sw, sy / sw]);
      },
    });
    assert(res.frames === 12 && progress === 12, `frames ${res.frames}, progress ${progress}`);
    assert(res.filename.endsWith(formatId === 'mp4-h264' ? '.mp4' : '.webm'), res.filename);
    for (const [cx, cy] of centroids) { near(cx, 160, 3, 'centroid x'); near(cy, 120, 3, 'centroid y'); }
    const url = URL.createObjectURL(res.blob);
    const video = document.createElement('video');
    video.muted = true; video.src = url;
    await videoEvent(video, 'loadeddata');
    assert(video.videoWidth === 320 && video.videoHeight === 240, `${video.videoWidth}x${video.videoHeight}`);
    near(video.duration, 12 / 24, 1.5 / 24, 'duration');
    URL.revokeObjectURL(url);
  });
}

await atest('record: overlay composited, cancel via AbortSignal', async () => {
  const state = recState({ overlay: true, spin: null });
  const ac = new AbortController();
  let frames = 0, lit = 0;
  const p = record({
    tables: TABLES, state, views: [{ element: state.element, visible: state.visible }], view: VIEW,
    width: 320, height: 240, fps: 30, duration: { seconds: 10 }, grain: 1, formatId: 'mp4-h264',
    onFrame: (canvas, i) => {
      if (i === 0) {   // tile label top-left: some bright pixels there
        const c = document.createElement('canvas'); c.width = 80; c.height = 80;
        const ctx = c.getContext('2d'); ctx.drawImage(canvas, 0, 0);
        const px = ctx.getImageData(0, 0, 80, 80).data;
        for (let j = 0; j < px.length; j += 4) if (px[j] + px[j + 1] + px[j + 2] > 300) lit++;
      }
    },
    onProgress: i => { frames = i; if (i === 5) ac.abort(); },
    signal: ac.signal,
  });
  let err = null;
  try { await p; } catch (e) { err = e; }
  assert(err && err.name === 'AbortError', `expected AbortError, got ${err}`);
  assert(frames === 5, `stopped after ${frames} frames`);
  assert(lit > 200, `overlay tile not drawn (${lit} bright px)`);
});

// ------------------------------------------------------------ report
const pass = results.filter(r => r.ok).length;
const out = document.getElementById('results');
for (const r of results) {
  const li = document.createElement('li');
  li.className = r.ok ? 'ok' : 'fail';
  li.textContent = (r.ok ? '✔ ' : '✘ ') + r.name + (r.ok ? '' : ' — ' + r.msg);
  out.append(li);
}
const summary = `${pass === results.length ? 'PASS' : 'FAIL'} ${pass}/${results.length}`;
document.getElementById('summary').textContent = summary;
document.title = summary;
console.log(`TESTS ${summary}` + results.filter(r => !r.ok).map(r => ` | FAIL ${r.name}: ${r.msg}`).join(''));
