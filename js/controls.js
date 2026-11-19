// Pointer-drag rotation, wheel zoom and two-finger pinch zoom.

import { qAxisAngle, qMul, qNormalize } from './math.js';

const RAD_PER_PX = Math.PI / 400;
export const ZOOM_MIN = 0.25;
export const ZOOM_MAX = 60;

export function clampZoom(z) {
  return Math.min(ZOOM_MAX, Math.max(ZOOM_MIN, z));
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
export function attachControls(el, { getOrientation, setOrientation, getZoom, setZoom }) {
  const pointers = new Map();     // pointerId -> [x, y]
  let pinch = null;               // { d0, z0 } while two pointers are down

  const distance = () => {
    const [a, b] = [...pointers.values()];
    return Math.hypot(a[0] - b[0], a[1] - b[1]);
  };

  el.addEventListener('pointerdown', e => {
    if (e.pointerType === 'mouse' && e.button !== 0) return;
    pointers.set(e.pointerId, [e.clientX, e.clientY]);
    el.setPointerCapture(e.pointerId);
    el.classList.add('dragging');
    if (pointers.size === 2) pinch = { d0: distance(), z0: getZoom() };
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
