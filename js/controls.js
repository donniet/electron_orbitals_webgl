// Pointer-drag rotation, wheel zoom and two-finger pinch zoom.

import { qAxisAngle, qMul, qNormalize } from './math.js';

const RAD_PER_PX = Math.PI / 400;
export const ZOOM_MIN = 0.25;
export const ZOOM_MAX = 200;

export function clampZoom(z) {
  return Math.min(ZOOM_MAX, Math.max(ZOOM_MIN, z));
}

// Fling: a drag released while still moving faster than this keeps spinning.
export const FLING_MIN_SPEED = 0.25;        // px / ms (~250 px/s)
export const SPIN_MAX_RATE = 0.006;         // rad / ms (~340 deg/s)
const VELOCITY_WINDOW = 80;                 // ms of recent motion used for the release velocity

// Spin {axis (view space), rate (rad/ms)} for a release velocity (vx, vy)
// in px/ms (screen y down), or null below the threshold. Same mapping as a
// drag: screen dx rotates about view +y, screen dy about view +x.
export function flingSpin(vx, vy, minSpeed = FLING_MIN_SPEED) {
  const speed = Math.hypot(vx, vy);
  if (!(speed >= minSpeed)) return null;
  return { axis: [vy / speed, vx / speed, 0], rate: Math.min(SPIN_MAX_RATE, speed * RAD_PER_PX) };
}

// Zoom after a pinch whose finger distance went from d0 to d1.
export function pinchZoom(startZoom, d0, d1) {
  if (!(d0 > 0) || !(d1 > 0)) return clampZoom(startZoom);
  return clampZoom(startZoom * (d1 / d0));
}

/**
 * Attach to `el`.
 *  - one pointer: rotate the atom in view space about the axis
 *    perpendicular to the drag (screen y is up in view space)
 *  - two pointers: pinch to zoom
 *  - mouse wheel: zoom
 */
export function attachControls(el, { getOrientation, setOrientation, getZoom, setZoom, setSpin }) {
  const pointers = new Map();     // pointerId -> [x, y]
  let pinch = null;               // { d0, z0 } while two pointers are down
  let pinched = false;            // a pinch happened during this gesture: no fling
  let samples = [];               // recent [t, x, y] of the rotating pointer

  const distance = () => {
    const [a, b] = [...pointers.values()];
    return Math.hypot(a[0] - b[0], a[1] - b[1]);
  };

  el.addEventListener('pointerdown', e => {
    if (e.pointerType === 'mouse' && e.button !== 0) return;
    pointers.set(e.pointerId, [e.clientX, e.clientY]);
    el.setPointerCapture(e.pointerId);
    el.classList.add('dragging');
    setSpin?.(null);                          // grabbing the atom stops any spin
    if (pointers.size === 1) { pinched = false; samples = [[e.timeStamp, e.clientX, e.clientY]]; }
    if (pointers.size === 2) { pinch = { d0: distance(), z0: getZoom() }; pinched = true; }
  });

  el.addEventListener('pointermove', e => {
    const last = pointers.get(e.pointerId);
    if (!last) return;
    const cur = [e.clientX, e.clientY];
    pointers.set(e.pointerId, cur);

    if (pointers.size >= 2) {
      if (pinch) setZoom(pinchZoom(pinch.z0, pinch.d0, distance()));
      return;
    }
    samples.push([e.timeStamp, cur[0], cur[1]]);
    while (samples.length > 2 && e.timeStamp - samples[0][0] > VELOCITY_WINDOW) samples.shift();
    const dx = cur[0] - last[0];
    const dy = cur[1] - last[1];
    const len = Math.hypot(dx, dy);
    if (len === 0) return;
    // Screen dx -> rotate about view +y; screen dy (down) -> about view +x.
    const dq = qAxisAngle([dy, dx, 0], len * RAD_PER_PX);
    setOrientation(qNormalize(qMul(dq, getOrientation())));
  });

  const end = e => {
    if (!pointers.has(e.pointerId)) return;
    if (pointers.size === 1 && !pinched && e.type === 'pointerup' && samples.length >= 2) {
      // Release velocity over the recent window; a pause before letting go
      // (no motion in the window) gives no fling.
      const [t0, x0, y0] = samples[0];
      const [t1, x1, y1] = samples[samples.length - 1];
      const dt = Math.max(e.timeStamp - t0, t1 - t0);
      if (dt > 0 && e.timeStamp - t1 < VELOCITY_WINDOW) setSpin?.(flingSpin((x1 - x0) / dt, (y1 - y0) / dt));
    }
    pointers.delete(e.pointerId);
    if (pointers.size < 2) pinch = null;   // remaining finger continues rotating from where it is
    if (pointers.size === 0) el.classList.remove('dragging');
    if (el.hasPointerCapture(e.pointerId)) el.releasePointerCapture(e.pointerId);
  };
  el.addEventListener('pointerup', end);
  el.addEventListener('pointercancel', end);

  el.addEventListener('wheel', e => {
    e.preventDefault();
    setZoom(clampZoom(getZoom() * Math.exp(-e.deltaY * 0.0015)));
  }, { passive: false });
}
