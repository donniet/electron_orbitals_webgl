// Element data and electron configurations.
//
// Effective nuclear charges are the Clementi & Raimondi (1963) SCF values.
// Orbitals are hydrogen-like with these Z_eff, in a fixed orientation in the
// atom frame. Within a p or d subshell orbitals fill by Hund's rule in a
// fixed order (p: x, y, z; d: xy, yz, xz, x²−y², z²): singly, then paired.

import { ORBITALS, SUBSHELLS, AUFBAU, orbitalsOf } from './orbitals.js';

export { ORBITALS } from './orbitals.js';

const CAPACITY = { 0: 2, 1: 6, 2: 10 };
const SUP = '⁰¹²³⁴⁵⁶⁷⁸⁹';
const sup = n => String(n).split('').map(c => SUP[+c]).join('');

// Ground-state subshell occupancy by the Aufbau order, with explicit
// exceptions (e.g. Cr, Cu) supplied per element.
export function aufbau(z) {
  const out = {};
  let left = z;
  for (const id of AUFBAU) {
    if (left <= 0) break;
    const cap = CAPACITY[SUBSHELLS.find(s => s.id === id).l];
    out[id] = Math.min(cap, left);
    left -= out[id];
  }
  return out;
}

// Distribute `count` electrons over the orbitals of a subshell (Hund).
export function hundFill(count, nOrbitals) {
  const occ = new Array(nOrbitals).fill(0);
  for (let e = 0; e < count; e++) occ[e % nOrbitals]++;
  return occ;
}

export function orbitalOccupancy(subshellOcc) {
  const occ = new Array(ORBITALS.length).fill(0);
  for (const [id, count] of Object.entries(subshellOcc)) {
    const idx = orbitalsOf(id);
    hundFill(count, idx.length).forEach((c, j) => { occ[idx[j]] = c; });
  }
  return occ;
}

export function electronConfig(subshellOcc) {
  return SUBSHELLS.filter(s => subshellOcc[s.id])
    .map(s => s.id + sup(subshellOcc[s.id])).join(' ');
}

const RAW = [
  // Z, sym, name, mass, group, category, electronegativity, 1st IE (eV), Z_eff
  [1,  'H',  'Hydrogen',  1.008,   1,  'Reactive nonmetal',     2.20, 13.598, { '1s': 1.000 }],
  [2,  'He', 'Helium',    4.0026,  18, 'Noble gas',             null, 24.587, { '1s': 1.688 }],
  [3,  'Li', 'Lithium',   6.94,    1,  'Alkali metal',          0.98, 5.392,  { '1s': 2.691, '2s': 1.279 }],
  [4,  'Be', 'Beryllium', 9.0122,  2,  'Alkaline earth metal',  1.57, 9.323,  { '1s': 3.685, '2s': 1.912 }],
  [5,  'B',  'Boron',     10.81,   13, 'Metalloid',             2.04, 8.298,  { '1s': 4.680, '2s': 2.576, '2p': 2.421 }],
  [6,  'C',  'Carbon',    12.011,  14, 'Reactive nonmetal',     2.55, 11.260, { '1s': 5.673, '2s': 3.217, '2p': 3.136 }],
  [7,  'N',  'Nitrogen',  14.007,  15, 'Reactive nonmetal',     3.04, 14.534, { '1s': 6.665, '2s': 3.847, '2p': 3.834 }],
  [8,  'O',  'Oxygen',    15.999,  16, 'Reactive nonmetal',     3.44, 13.618, { '1s': 7.658, '2s': 4.492, '2p': 4.453 }],
  [9,  'F',  'Fluorine',  18.998,  17, 'Halogen',               3.98, 17.423, { '1s': 8.650, '2s': 5.128, '2p': 5.100 }],
  [10, 'Ne', 'Neon',      20.180,  18, 'Noble gas',             null, 21.565, { '1s': 9.642, '2s': 5.758, '2p': 5.758 }],
];

function periodOf(z) {
  return z <= 2 ? 1 : z <= 10 ? 2 : z <= 18 ? 3 : 4;
}

export const ELEMENTS = RAW.map(([z, symbol, name, mass, group, category, en, ie, zeff, configOverride]) => {
  const subshells = configOverride || aufbau(z);
  return {
    z, symbol, name, mass, group, category,
    period: periodOf(z),
    electronegativity: en,
    ionizationEnergy: ie,
    subshells,                                         // { '1s': 2, ... }
    occupancy: orbitalOccupancy(subshells),            // per ORBITALS index
    zeff: ORBITALS.map(o => zeff[o.subshell] ?? 0),    // per ORBITALS index
    zeffBySubshell: zeff,
    config: electronConfig(subshells),
  };
});

export function elementByZ(z) {
  return ELEMENTS[z - 1];
}
