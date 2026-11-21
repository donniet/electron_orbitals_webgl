// Video recording: renders the current view offline, frame by frame, at a
// chosen resolution and frame rate, encodes each frame as it is rendered
// (WebCodecs VideoEncoder) and writes an MP4 or WebM file.
//
// Frames are never stored raw (4K at 60 fps is ~2 GB per second); only the
// encoded chunks are kept until the file is written.
//
// The output frame is a centre crop of the on-screen view: the largest
// rectangle of the output's aspect ratio that fits in the view, rendered at
// the output resolution (same bohr-per-pixel scale as the crop, so what is
// inside the crop rectangle is what the video shows). With two atoms the
// output is split into its own panes, each centred on its atom (paneCrops
// gives the matching rectangle in each on-screen pane).

import { createRenderer, splitViews } from './renderer.js';
import { qAxisAngle, qMul, qNormalize } from './math.js';
import { drawOverlay } from './ui.js';
import { muxMp4 } from './mux-mp4.js';
import { muxWebm } from './mux-webm.js';

// ------------------------------------------------------------ options
// Output sizes offered besides "current" (the view's own device-pixel size).
export const RESOLUTIONS = [
  { id: '1280x720', w: 1280, h: 720, label: '1280 × 720 (720p)' },
  { id: '1920x1080', w: 1920, h: 1080, label: '1920 × 1080 (1080p)' },
  { id: '2560x1440', w: 2560, h: 1440, label: '2560 × 1440 (1440p)' },
  { id: '3840x2160', w: 3840, h: 2160, label: '3840 × 2160 (4K)' },
  { id: '1080x1080', w: 1080, h: 1080, label: '1080 × 1080 (square)' },
  { id: '1080x1920', w: 1080, h: 1920, label: '1080 × 1920 (portrait)' },
];

export const FPS_CHOICES = [24, 30, 50, 60];

// Bits per output pixel per frame. Fresh noise in every frame defeats the
// codecs' motion prediction, so these are high by usual video standards.
// Calibrated by eye and PSNR on source vs decoded frames (1280x720, 60 fps):
// at 0.3 the dense core smears into smooth blobs; 0.8 with 2 px grain is the
// first setting that keeps the grain and its colours.
export const QUALITY = {
  standard: { label: 'Standard', bpp: 0.3 },
  high: { label: 'High', bpp: 0.8 },
  max: { label: 'Max', bpp: 1.6 },
};

// Keyframe interval (seconds): seeking granularity.
const KEY_INTERVAL_S = 2;

// The view's size rounded down to even numbers (4:2:0 video needs even sizes).
export function evenSize(w, h) {
  return [Math.max(2, w & ~1), Math.max(2, h & ~1)];
}

// Largest centred rectangle with the output's aspect ratio inside a W x H
// view (device px, top-left origin).
export function cropRect(W, H, outW, outH) {
  const s = Math.min(W / outW, H / outH);
  const w = outW * s, h = outH * s;
  return { x: (W - w) / 2, y: (H - h) / 2, w, h, scale: s };   // scale: view px per output px
}

// What each output pane shows, as rectangles in the on-screen view (device
// px, top-left origin). Output pane j is centred on its atom like on-screen
// pane j, at the recording's scale (view px per output px), so its rectangle
// is centred on that pane with the output pane's size times the scale.
// panes: the on-screen pane rects (GL origin, bottom-left), as renderViews uses.
export function paneCrops(W, H, panes, outW, outH, grain = 1) {
  const { scale } = cropRect(W, H, outW, outH);
  return splitViews(outW, outH, grain, panes.length).map((r, j) => {
    const p = panes[j];
    const cx = p.x + p.w / 2, cy = H - (p.y + p.h / 2);
    const w = r.w * scale, h = r.h * scale;
    return { x: cx - w / 2, y: cy - h / 2, w, h };
  });
}

// ------------------------------------------------------------ frame plan
// spin: { axis, rate (rad/ms) } in view space, or null.
// Rotation period (seconds), or null without spin.
export function rotationPeriod(spin) {
  return spin && spin.rate > 0 ? (2 * Math.PI) / spin.rate / 1000 : null;
}

/**
 * Frames to record. duration: { rotation: true } for exactly one turn, or
 * { seconds }. One turn uses n = round(period * fps) frames and a step of
 * exactly 2 pi / n, so the video loops seamlessly (frame n would be frame 0).
 * Returns { n, step (rad per frame), axis, orientation(i) }.
 */
export function framePlan(orientation, spin, fps, duration) {
  let n, step = 0;
  const axis = spin ? spin.axis : [0, 1, 0];
  if (duration.rotation) {
    const period = rotationPeriod(spin);
    if (!period) throw new Error('One rotation needs a spinning atom.');
    n = Math.max(1, Math.round(period * fps));
    step = (2 * Math.PI) / n;
  } else {
    n = Math.max(1, Math.round(duration.seconds * fps));
    step = spin ? (spin.rate * 1000) / fps : 0;
  }
  // Same composition as the live spin in main.js: rotate in view space.
  const at = i => qNormalize(qMul(qAxisAngle(axis, step * i), orientation));
  return { n, step, axis, orientation: at };
}

// ------------------------------------------------------------ codecs
// Levels: the lowest that allows the frame size, sample rate and bitrate,
// so the stream declares limits it keeps (hardware decoders check them).
// If none fits, the last; the encoder probe then decides.

// H.264 High profile (ITU-T H.264 Table A-1; High MaxBR = 1.25 x the table).
const H264_LEVELS = [   // [level_idc, max MB per second, max frame MBs, max kbit/s]
  [30, 40500, 1620, 12500], [31, 108000, 3600, 17500], [32, 216000, 5120, 25000],
  [40, 245760, 8192, 25000], [41, 245760, 8192, 62500], [42, 522240, 8704, 62500],
  [50, 589824, 22080, 168750], [51, 983040, 36864, 300000], [52, 2073600, 36864, 300000],
  [60, 4177920, 139264, 300000], [61, 8355840, 139264, 600000], [62, 16711680, 139264, 1000000],
];
export function h264Codec(w, h, fps, bitrate = 0) {
  const fs = Math.ceil(w / 16) * Math.ceil(h / 16);
  const lv = H264_LEVELS.find(([, mbps, mfs, br]) => fs <= mfs && fs * fps <= mbps && bitrate <= br * 1000)
    ?? H264_LEVELS.at(-1);
  return `avc1.6400${lv[0].toString(16).toUpperCase().padStart(2, '0')}`;
}

// VP9 profile 0 (WebM project level table).
const VP9_LEVELS = [   // [level, max luma samples / s, max picture size, max Mbit/s]
  [31, 36864000, 983040, 12], [40, 83558400, 2228224, 18], [41, 160432128, 2228224, 30],
  [50, 311951360, 8912896, 60], [51, 588251136, 8912896, 120], [52, 1176502272, 8912896, 180],
  [60, 1176502272, 35651584, 180], [61, 2353004544, 35651584, 240], [62, 4706009088, 35651584, 480],
];
export function vp9Codec(w, h, fps, bitrate = 0) {
  const lv = VP9_LEVELS.find(([, rate, size, br]) => w * h <= size && w * h * fps <= rate && bitrate <= br * 1e6)
    ?? VP9_LEVELS.at(-1);
  return `vp09.00.${lv[0]}.08`;
}

// AV1 main profile; main tier (M), or high tier (H) where the bitrate needs it.
const AV1_LEVELS = [   // [seq_level_idx, max display rate, max picture size, main Mbit/s, high Mbit/s]
  [5, 24969600, 665856, 10, 0], [8, 66846720, 2359296, 12, 30], [9, 133693440, 2359296, 20, 50],
  [12, 267386880, 8912896, 30, 100], [13, 534773760, 8912896, 40, 160], [14, 1069547520, 8912896, 60, 240],
  [16, 1069547520, 35651584, 60, 240], [17, 2139095040, 35651584, 100, 480], [18, 4278190080, 35651584, 160, 800],
];
export function av1Codec(w, h, fps, bitrate = 0) {
  for (const [idx, rate, size, main, high] of AV1_LEVELS) {
    if (w * h > size || w * h * fps > rate) continue;
    if (bitrate <= main * 1e6) return `av01.0.${String(idx).padStart(2, '0')}M.08`;
    if (bitrate <= high * 1e6) return `av01.0.${String(idx).padStart(2, '0')}H.08`;
  }
  return 'av01.0.18H.08';
}

// Formats in order of preference (most widely playable first).
export const FORMATS = [
  { id: 'mp4-h264', label: 'MP4 · H.264', mime: 'video/mp4', ext: 'mp4', tag: 'h264', codec: h264Codec },
  { id: 'webm-vp9', label: 'WebM · VP9', mime: 'video/webm', ext: 'webm', tag: 'vp9', codec: vp9Codec },
  { id: 'webm-av1', label: 'WebM · AV1', mime: 'video/webm', ext: 'webm', tag: 'av1', codec: av1Codec },
  { id: 'webm-vp8', label: 'WebM · VP8', mime: 'video/webm', ext: 'webm', tag: 'vp8', codec: () => 'vp8' },
];

export function bitrateFor(w, h, fps, quality) {
  return Math.round(QUALITY[quality].bpp * w * h * fps);
}

export function estimateBytes(w, h, fps, n, quality) {
  return (bitrateFor(w, h, fps, quality) * (n / fps)) / 8;
}

function encoderConfig(format, w, h, fps, quality) {
  const config = {
    codec: format.codec(w, h, fps, bitrateFor(w, h, fps, quality)), width: w, height: h, framerate: fps,
    bitrate: bitrateFor(w, h, fps, quality), latencyMode: 'quality',
  };
  if (format.id === 'mp4-h264') config.avc = { format: 'avc' };
  return config;
}

export function videoSupported() {
  return typeof VideoEncoder === 'function' && typeof VideoFrame === 'function';
}

// The formats this browser can encode at these settings.
export async function supportedFormats(w, h, fps, quality) {
  if (!videoSupported()) return [];
  const ok = await Promise.all(FORMATS.map(async f => {
    try { return (await VideoEncoder.isConfigSupported(encoderConfig(f, w, h, fps, quality))).supported; }
    catch { return false; }
  }));
  return FORMATS.filter((f, i) => ok[i]);
}

// ------------------------------------------------------------ recording
/**
 * Render and encode a video. Options:
 *   tables, state     the app's tables and (live) state; copied, not changed
 *   views             [{ element, visible }] panes, as in main.js
 *   view              { width, height, bohrPerPixel, dpr }: the on-screen view
 *                     (device px), its scale and device pixel ratio
 *   width, height     output size (even)
 *   fps, duration     see framePlan
 *   grain             1 or 2: output pixels per noise sample (render scale k)
 *   formatId, quality see FORMATS, QUALITY
 *   onProgress(i, n)  after each frame
 *   onFrame(canvas, i) optional: sees each source frame before it is encoded
 *   signal            AbortSignal to cancel
 * Resolves to { blob, filename, frames, bytes }.
 */
export async function record(opts) {
  const { tables, state, views, view, width, height, fps, duration, grain = 2,
    formatId, quality = 'high', onProgress = () => {}, onFrame = null, signal } = opts;
  const format = FORMATS.find(f => f.id === formatId);
  if (!format) throw new Error(`Unknown format ${formatId}`);
  const plan = framePlan(state.orientation, state.spin, fps, duration);
  const crop = cropRect(view.width, view.height, width, height);

  const glCanvas = document.createElement('canvas');
  glCanvas.width = width; glCanvas.height = height;
  const renderer = createRenderer(glCanvas, tables);
  let source = glCanvas, compose = null, overlay = null;
  if (state.overlay) {
    overlay = document.createElement('canvas');
    compose = document.createElement('canvas');
    overlay.width = compose.width = width;
    overlay.height = compose.height = height;
    source = compose;
  }

  const chunks = [];
  let description = null, failure = null;
  const encoder = new VideoEncoder({
    output: (chunk, meta) => {
      const data = new Uint8Array(chunk.byteLength);
      chunk.copyTo(data);
      chunks.push({ data, pts: chunk.timestamp, key: chunk.type === 'key' });
      if (meta?.decoderConfig?.description) description = meta.decoderConfig.description;
    },
    error: e => { failure = e; },
  });
  const config = encoderConfig(format, width, height, fps, quality);
  encoder.configure(config);

  const keyEvery = Math.max(1, Math.round(KEY_INTERVAL_S * fps));
  const frameUs = 1e6 / fps;
  const rects = splitViews(width, height, grain, views.length);
  const rs = {
    ...state,
    bohrPerPixel: view.bohrPerPixel * crop.scale,
    baseHalf: null,
  };
  const overlayScale = view.dpr / crop.scale;   // keep the overlay's on-screen proportions
  try {
    for (let i = 0; i < plan.n; i++) {
      if (signal?.aborted) throw new DOMException('Recording cancelled', 'AbortError');
      if (failure) throw failure;
      rs.orientation = plan.orientation(i);
      rs.frame = i;
      rs.seed = state.debugSeed ? i + 1 : (Math.random() * 4294967296) >>> 0;
      const geos = renderer.renderViews(rs, grain, views.map((v, j) => ({ ...v, rect: rects[j] })));
      if (compose) {
        drawOverlay(overlay, rs, views.map((v, j) => ({ element: v.element, geo: geos[j] })), overlayScale);
        const ctx = compose.getContext('2d');
        ctx.drawImage(glCanvas, 0, 0);
        ctx.drawImage(overlay, 0, 0);
      }
      if (onFrame) onFrame(source, i);
      const frame = new VideoFrame(source, { timestamp: Math.round(i * frameUs), duration: Math.round(frameUs) });
      encoder.encode(frame, { keyFrame: i % keyEvery === 0 });
      frame.close();
      // Backpressure: keep only a few frames queued in the encoder.
      while (encoder.encodeQueueSize > 2 && !failure) {
        await new Promise(r => {                           // (timeout: never hang on an encoder error)
          encoder.addEventListener('dequeue', r, { once: true });
          setTimeout(r, 50);
        });
      }
      onProgress(i + 1, plan.n);
      await new Promise(r => setTimeout(r));               // let the page update
    }
    await encoder.flush();
    if (failure) throw failure;
  } finally {
    if (encoder.state !== 'closed') encoder.close();
    renderer.gl.getExtension('WEBGL_lose_context')?.loseContext();
  }

  let bytes;
  if (format.mime === 'video/mp4') {
    if (!description) throw new Error('The H.264 encoder gave no avcC description.');
    const timescale = fps * 1000;
    bytes = muxMp4({ width, height, timescale, frameDuration: 1000, description,
      chunks: chunks.map(c => ({ ...c, pts: Math.round((c.pts * timescale) / 1e6) })) });
  } else {
    bytes = muxWebm({ width, height, fps, codec: config.codec, description, chunks });
  }
  const names = views.map(v => v.element.symbol).join('-');
  const filename = `orbitals-${names}-${width}x${height}-${fps}fps-${format.tag}.${format.ext}`;
  return { blob: new Blob([bytes], { type: format.mime }), filename, frames: plan.n, bytes: bytes.length };
}
