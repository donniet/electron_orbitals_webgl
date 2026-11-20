// The set of selected atoms: A (the primary) and an optional B shown side by
// side. Everything that should cover "both atoms" (the orbital legend, the
// Valence filter, which shells start expanded, which periodic-table rows the
// collapsed strip shows) is derived here from the union of the selection, so
// the single-atom case is just a selection of one.
//
// Orbital visibility stays one array per orbital (state.visible): an orbital
// slot such as 2px is the same for every element, so hiding it hides it in
// both panes.

import { ORBITALS, isValence, valenceSubshells } from './orbitals.js';

/** [A] or [A, B]. */
export function atomsOf(state) {
  return state.second ? [state.element, state.second] : [state.element];
}

/** Subshell ids occupied in any of the atoms. */
export function subshellUnion(atoms) {
  const ids = new Set();
  for (const el of atoms) for (const id of Object.keys(el.subshells)) ids.add(id);
  return ids;
}

/** Is orbital i occupied in any of the atoms? */
export function occupiedAny(atoms, i) {
  return atoms.some(el => el.occupancy[i] > 0);
}

/** Orbital indices that are valence in any of the atoms. */
export function valenceUnion(atoms) {
  const out = new Set();
  ORBITALS.forEach((o, i) => { if (atoms.some(el => isValence(o, el))) out.add(i); });
  return out;
}

/** Highest period among the atoms (the legend lists subshells up to it). */
export function maxPeriod(atoms) {
  return Math.max(...atoms.map(el => el.period));
}

/**
 * Shells (n) expanded by default in the legend: every shell when the atoms
 * have few subshells between them, otherwise the shells holding valence
 * subshells of any atom.
 */
export function defaultExpandedShells(atoms) {
  const occupied = subshellUnion(atoms);
  if (occupied.size <= 6) return new Set([...occupied].map(id => +id[0]));
  const out = new Set();
  for (const el of atoms) for (const id of valenceSubshells(el)) out.add(+id[0]);
  return out;
}

/**
 * Rows shown by the collapsed periodic table: each row in `fullRows` in
 * full, plus thin "peek" rows for the neighbours of a full row (in table
 * order). Returns [{ row, peek, gridRow }] in display order.
 */
export function collapsedRows(allRows, fullRows) {
  const rows = [...new Set(allRows)].sort((a, b) => a - b);
  const full = new Set(fullRows);
  const shown = new Map();                         // row -> peek?
  rows.forEach((r, i) => {
    if (!full.has(r)) return;
    shown.set(r, false);
    for (const j of [i - 1, i + 1]) {
      if (j >= 0 && j < rows.length && !full.has(rows[j]) && !shown.has(rows[j])) shown.set(rows[j], true);
    }
  });
  return [...shown.entries()].sort((a, b) => a[0] - b[0])
    .map(([row, peek], i) => ({ row, peek, gridRow: i + 1 }));
}

/** Periodic-table row of an element (the f-block has its own row, 8). */
export function tableRow(el) {
  return el.group == null ? 8 : el.period;
}

