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

const AR_CORE = { '1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6 };
function core(s1, s2, p2, s3, p3) {
  return { '1s': s1, '2s': s2, '2p': p2, '3s': s3, '3p': p3 };
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
  // Period 4. Z_eff from the same table. Cr and Cu take one 4s electron
  // into 3d (half-filled / filled d), overriding the Aufbau order.
  [19, 'K',  'Potassium',  39.098, 1,  'Alkali metal',          0.82, 4.341,  { ...core(18.490, 13.006, 15.027, 8.680, 7.726), '4s': 3.495 }],
  [20, 'Ca', 'Calcium',    40.078, 2,  'Alkaline earth metal',  1.00, 6.113,  { ...core(19.473, 13.776, 16.041, 9.602, 8.658), '4s': 4.398 }],
  [21, 'Sc', 'Scandium',   44.956, 3,  'Transition metal',      1.36, 6.561,  { ...core(20.457, 14.574, 17.055, 10.340, 9.406), '3d': 7.120, '4s': 4.632 }],
  [22, 'Ti', 'Titanium',   47.867, 4,  'Transition metal',      1.54, 6.828,  { ...core(21.441, 15.377, 18.065, 11.033, 10.104), '3d': 8.141, '4s': 4.817 }],
  [23, 'V',  'Vanadium',   50.942, 5,  'Transition metal',      1.63, 6.746,  { ...core(22.426, 16.181, 19.073, 11.709, 10.785), '3d': 8.983, '4s': 4.981 }],
  [24, 'Cr', 'Chromium',   51.996, 6,  'Transition metal',      1.66, 6.767,  { ...core(23.414, 16.984, 20.075, 12.368, 11.466), '3d': 9.757, '4s': 5.133 },
    { ...AR_CORE, '3d': 5, '4s': 1 }],
  [25, 'Mn', 'Manganese',  54.938, 7,  'Transition metal',      1.55, 7.434,  { ...core(24.396, 17.794, 21.084, 13.018, 12.109), '3d': 10.528, '4s': 5.283 }],
  [26, 'Fe', 'Iron',       55.845, 8,  'Transition metal',      1.83, 7.902,  { ...core(25.381, 18.599, 22.089, 13.676, 12.778), '3d': 11.180, '4s': 5.434 }],
  [27, 'Co', 'Cobalt',     58.933, 9,  'Transition metal',      1.88, 7.881,  { ...core(26.367, 19.405, 23.092, 14.322, 13.435), '3d': 11.855, '4s': 5.576 }],
  [28, 'Ni', 'Nickel',     58.693, 10, 'Transition metal',      1.91, 7.640,  { ...core(27.353, 20.213, 24.095, 14.961, 14.085), '3d': 12.530, '4s': 5.711 }],
  [29, 'Cu', 'Copper',     63.546, 11, 'Transition metal',      1.90, 7.726,  { ...core(28.339, 21.020, 25.097, 15.594, 14.731), '3d': 13.201, '4s': 5.842 },
    { ...AR_CORE, '3d': 10, '4s': 1 }],
  [30, 'Zn', 'Zinc',       65.38,  12, 'Transition metal',      1.65, 9.394,  { ...core(29.325, 21.828, 26.098, 16.219, 15.369), '3d': 13.878, '4s': 5.965 }],
  [31, 'Ga', 'Gallium',    69.723, 13, 'Post-transition metal', 1.81, 5.999,  { ...core(30.309, 22.599, 27.091, 16.996, 16.204), '3d': 15.093, '4s': 7.067, '4p': 6.222 }],
  [32, 'Ge', 'Germanium',  72.630, 14, 'Metalloid',             2.01, 7.900,  { ...core(31.294, 23.365, 28.082, 17.790, 17.014), '3d': 16.251, '4s': 8.044, '4p': 6.780 }],
  [33, 'As', 'Arsenic',    74.922, 15, 'Metalloid',             2.18, 9.789,  { ...core(32.278, 24.127, 29.074, 18.596, 17.850), '3d': 17.378, '4s': 8.944, '4p': 7.449 }],
  [34, 'Se', 'Selenium',   78.971, 16, 'Reactive nonmetal',     2.55, 9.752,  { ...core(33.262, 24.888, 30.065, 19.403, 18.705), '3d': 18.477, '4s': 9.758, '4p': 8.287 }],
  [35, 'Br', 'Bromine',    79.904, 17, 'Halogen',               2.96, 11.814, { ...core(34.247, 25.643, 31.056, 20.219, 19.571), '3d': 19.559, '4s': 10.553, '4p': 9.028 }],
  [36, 'Kr', 'Krypton',    83.798, 18, 'Noble gas',             3.00, 14.000, { ...core(35.232, 26.398, 32.047, 21.033, 20.434), '3d': 20.626, '4s': 11.316, '4p': 9.338 }],
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
