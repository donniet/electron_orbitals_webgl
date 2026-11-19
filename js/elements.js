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

// Noble-gas cores for shorthand configurations (largest first).
const NOBLE_CORES = [
  ['Ar', 18, { '1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6 }],
  ['Ne', 10, { '1s': 2, '2s': 2, '2p': 6 }],
  ['He', 2, { '1s': 2 }],
];

// "[Ne] 3s² 3p⁴" for Z > 10; the full configuration for periods 1-2,
// which are short enough to read in full.
export function shortConfig(z, subshellOcc) {
  if (z <= 10) return electronConfig(subshellOcc);
  for (const [sym, zc, core] of NOBLE_CORES) {
    if (z <= zc) continue;
    if (!Object.entries(core).every(([id, c]) => subshellOcc[id] === c)) continue;
    const rest = Object.fromEntries(Object.entries(subshellOcc).filter(([id]) => !(id in core)));
    return `[${sym}] ` + electronConfig(rest);
  }
  return electronConfig(subshellOcc);
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
  // Period 3. Z_eff: Clementi & Raimondi (1963), as tabulated on Wikipedia
  // "Effective nuclear charge" (fetched 2026-11-18).
  [11, 'Na', 'Sodium',     22.990, 1,  'Alkali metal',          0.93, 5.139,  { '1s': 10.626, '2s': 6.571, '2p': 6.802, '3s': 2.507 }],
  [12, 'Mg', 'Magnesium',  24.305, 2,  'Alkaline earth metal',  1.31, 7.646,  { '1s': 11.609, '2s': 7.392, '2p': 7.826, '3s': 3.308 }],
  [13, 'Al', 'Aluminium',  26.982, 13, 'Post-transition metal', 1.61, 5.986,  { '1s': 12.591, '2s': 8.214, '2p': 8.963, '3s': 4.117, '3p': 4.066 }],
  [14, 'Si', 'Silicon',    28.085, 14, 'Metalloid',             1.90, 8.152,  { '1s': 13.575, '2s': 9.020, '2p': 9.945, '3s': 4.903, '3p': 4.285 }],
  [15, 'P',  'Phosphorus', 30.974, 15, 'Reactive nonmetal',     2.19, 10.487, { '1s': 14.558, '2s': 9.825, '2p': 10.961, '3s': 5.642, '3p': 4.886 }],
  [16, 'S',  'Sulfur',     32.06,  16, 'Reactive nonmetal',     2.58, 10.360, { '1s': 15.541, '2s': 10.629, '2p': 11.977, '3s': 6.367, '3p': 5.482 }],
  [17, 'Cl', 'Chlorine',   35.45,  17, 'Halogen',               3.16, 12.968, { '1s': 16.524, '2s': 11.430, '2p': 12.993, '3s': 7.068, '3p': 6.116 }],
  [18, 'Ar', 'Argon',      39.948, 18, 'Noble gas',             null, 15.760, { '1s': 17.508, '2s': 12.230, '2p': 14.008, '3s': 7.757, '3p': 6.764 }],
];

export function periodOf(z) {
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
    configShort: shortConfig(z, subshells),
  };
});

export function elementByZ(z) {
  return ELEMENTS[z - 1];
}
