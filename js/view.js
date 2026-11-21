// View defaults and smooth view-scale changes (e.g. when switching elements).
//
// The animator eases a positive scale value toward a target in log space,
// so growing and shrinking by the same factor take the same visual pace.
// Retargeting mid-animation starts from the value currently on screen, so
// rapid changes chain without jumps.

import { qAxisAngle, qMul } from './math.js';

// Default orientation: tilted so x, y and z lobes are all distinguishable.
export function defaultOrientation() {
  return qMul(qAxisAngle([1, 0, 0], 0.45), qAxisAngle([0, 1, 0], -0.6));
}

const easeInOutCubic = t => (t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2);

export function prefersReducedMotion() {
  return typeof matchMedia === 'function' && matchMedia('(prefers-reduced-motion: reduce)').matches;
}

export function createScaleAnimator({ duration = 450 } = {}) {
  let from = null, to = null, current = null, t0 = 0;

  const anim = {
    animating: false,
    /** Current value for `target` at time `now` (ms). */
    value(target, now) {
      if (current === null || duration <= 0) {
        from = to = current = target;
        anim.animating = false;
        return current;
      }
      if (target !== to) {             // new target: start from what is shown
        from = current;
        to = target;
        t0 = now;
      }
      const t = Math.min(1, Math.max(0, (now - t0) / duration));
      current = t >= 1 ? to : Math.exp(Math.log(from) + (Math.log(to) - Math.log(from)) * easeInOutCubic(t));
      anim.animating = t < 1;
      return current;
    },
  };
  return anim;
}
