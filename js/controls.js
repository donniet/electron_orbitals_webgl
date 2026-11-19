// Pointer-drag rotation and wheel zoom.

import { qAxisAngle, qMul, qNormalize } from './math.js';

const RAD_PER_PX = Math.PI / 400;

/**
 * Attach to `el`. getState/setState operate on { orientation, zoom }.
 * Dragging rotates the atom in view space about the axis perpendicular to
 * the drag (screen y is up in view space, so drag-down tilts toward the
 * viewer from the top).
 */
export function attachControls(el, { getOrientation, setOrientation, getZoom, setZoom }) {
  let last = null;

  el.addEventListener('pointerdown', e => {
    if (e.button !== 0) return;
    last = [e.clientX, e.clientY];
    el.setPointerCapture(e.pointerId);
    el.classList.add('dragging');
  });

  el.addEventListener('pointermove', e => {
    if (!last) return;
    const dx = e.clientX - last[0];
    const dy = e.clientY - last[1];
    last = [e.clientX, e.clientY];
    const len = Math.hypot(dx, dy);
    if (len === 0) return;
    // Screen dx -> rotate about view +y; screen dy (down) -> about view +x.
    const dq = qAxisAngle([dy, dx, 0], len * RAD_PER_PX);
    setOrientation(qNormalize(qMul(dq, getOrientation())));
  });

  const end = e => {
    if (!last) return;
    last = null;
    el.classList.remove('dragging');
    if (el.hasPointerCapture(e.pointerId)) el.releasePointerCapture(e.pointerId);
  };
  el.addEventListener('pointerup', end);
  el.addEventListener('pointercancel', end);

  el.addEventListener('wheel', e => {
    e.preventDefault();
    const z = getZoom() * Math.exp(-e.deltaY * 0.0015);
    setZoom(Math.min(60, Math.max(0.25, z)));
  }, { passive: false });
}
