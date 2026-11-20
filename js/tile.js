// Periodic-table tile (category colour, Z, symbol, atomic mass), drawn with
// the 2D canvas API so the same rendering is used in the side-panel card
// and on the canvas overlay (and therefore in saved PNGs).
//
// Category colours stay defined once, as CSS custom properties
// (--cat-<slug> in css/style.css); this module reads them from there.

export const TILE_W = 84;   // CSS px at scale 1
export const TILE_H = 96;

export function categorySlug(category) {
  return category.toLowerCase().replace(/[^a-z]+/g, '-');
}

export function categoryClass(category) {
  return 'cat-' + categorySlug(category);
}

export function categoryColor(category) {
  const v = getComputedStyle(document.documentElement)
    .getPropertyValue('--cat-' + categorySlug(category)).trim();
  return v || '#555b66';
}

/** Draw the tile for `el` with its top-left at (x, y), in canvas pixels. */
export function drawElementTile(ctx, el, x, y, scale = 1) {
  const w = TILE_W * scale, h = TILE_H * scale;
  ctx.save();
  ctx.fillStyle = categoryColor(el.category);
  ctx.beginPath();
  ctx.roundRect(x, y, w, h, 6 * scale);
  ctx.fill();

  ctx.fillStyle = '#fff';
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';
  ctx.font = `${12 * scale}px system-ui, sans-serif`;
  ctx.fillText(String(el.z), x + 7 * scale, y + 5 * scale);

  ctx.textAlign = 'center';
  ctx.textBaseline = 'alphabetic';
  ctx.font = `700 ${38 * scale}px system-ui, sans-serif`;
  ctx.fillText(el.symbol, x + w / 2, y + 58 * scale);

  ctx.globalAlpha = 0.9;
  ctx.font = `${11 * scale}px system-ui, sans-serif`;
  ctx.fillText(String(Number(el.mass.toFixed(3))), x + w / 2, y + 78 * scale);
  ctx.restore();
}

/** Size `canvas` for the device pixel ratio and draw the tile into it. */
export function renderTileCanvas(canvas, el, scale = 1) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  canvas.width = Math.round(TILE_W * scale * dpr);
  canvas.height = Math.round(TILE_H * scale * dpr);
  canvas.style.width = `${TILE_W * scale}px`;
  canvas.style.height = `${TILE_H * scale}px`;
  const ctx = canvas.getContext('2d');
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  drawElementTile(ctx, el, 0, 0, scale * dpr);
}
