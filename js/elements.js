// Element data for H..Ne and the orbitals we render.
//
// Effective nuclear charges are the Clementi & Raimondi (1963) SCF values.
// Orbitals are hydrogen-like with these Z_eff, in a fixed orientation in the
// atom frame. p orbitals fill by Hund's rule: px, py, pz singly, then paired.

// Orbital render order (index into shader uniform arrays).
// kind: 's1' | 's2' | 'p' ; axis is the lobe axis for p orbitals.
// label is HTML; hue is in degrees.
export const ORBITALS = [
  { id: '1s',  label: '1s',  n: 1, l: 0, kind: 's1', axis: -1, hue: 0,   zeffKey: '1s' },
  { id: '2s',  label: '2s',  n: 2, l: 0, kind: 's2', axis: -1, hue: 45,  zeffKey: '2s' },
  { id: '2px', label: '2p<sub>x</sub>', n: 2, l: 1, kind: 'p',  axis: 0,  hue: 195, zeffKey: '2p' },
  { id: '2py', label: '2p<sub>y</sub>', n: 2, l: 1, kind: 'p',  axis: 1,  hue: 275, zeffKey: '2p' },
  { id: '2pz', label: '2p<sub>z</sub>', n: 2, l: 1, kind: 'p',  axis: 2,  hue: 125, zeffKey: '2p' },
];

const SUP = { 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶' };

const P_FILL = [ // Hund's rule occupancy of (px, py, pz) for 0..6 electrons
  [0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1], [2, 1, 1], [2, 2, 1], [2, 2, 2],
];

function occupancy(z) {
  const s1 = Math.min(z, 2);
  const s2 = Math.min(Math.max(z - 2, 0), 2);
  const p = Math.max(z - 4, 0);
  return [s1, s2, ...P_FILL[p]];
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

export const ELEMENTS = RAW.map(([z, symbol, name, mass, group, category, en, ie, zeff]) => {
  const occ = occupancy(z);
  return {
    z, symbol, name, mass, group, category,
    period: z <= 2 ? 1 : 2,
    electronegativity: en,
    ionizationEnergy: ie,
    occupancy: occ,                                   // per ORBITALS index
    zeff: ORBITALS.map(o => zeff[o.zeffKey] ?? 0),    // per ORBITALS index
    config: electronConfig(occ),
  };
});

function electronConfig(occ) {
  const parts = [];
  if (occ[0]) parts.push('1s' + SUP[occ[0]]);
  if (occ[1]) parts.push('2s' + SUP[occ[1]]);
  const p = occ[2] + occ[3] + occ[4];
  if (p) parts.push('2p' + SUP[p]);
  return parts.join(' ');
}

export function elementByZ(z) {
  return ELEMENTS[z - 1];
}
