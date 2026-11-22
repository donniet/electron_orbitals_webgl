// Loads the moment tables written by tools/gen_tables.py:
//   data/orbital_tables.json      metadata (N, uMax, rows, channels, odd)
//   data/orbital_tables.bin       float32 little-endian, [row][i][RGBA]
//   data/orbital_tables_odd.bin   same layout, odd half-line moments
//   data/orbital_tables_radial.bin  float32, [row][i], q2 = R^2 / u^(2l)
//   data/orbital_tables_momentum.bin  float32, [row][i], scaled momentum
//                                    q = p n / Z at cumulative F = (i + 0.5) / N
// Texels hold (I0, I2, I4, I6) and (J1, J3, J5, 0) of one radial function
// at u = i * uMax / (N - 1); rows[r].extent is the radius (reference units)
// beyond which row r's density is negligible.

const DEFAULT_BASE = new URL('../data/', import.meta.url);

// Always revalidate the data files (a 304 when unchanged): a cached old
// metadata file next to new tables or code would load mismatched data.
const FETCH_OPTS = { cache: 'no-cache' };

export async function loadTables(base = DEFAULT_BASE) {
  const [metaRes, binRes] = await Promise.all([
    fetch(new URL('orbital_tables.json', base), FETCH_OPTS),
    fetch(new URL('orbital_tables.bin', base), FETCH_OPTS),
  ]);
  if (!metaRes.ok || !binRes.ok) {
    throw new Error(`Could not load orbital tables (${metaRes.status}/${binRes.status}). ` +
      'Run python3 tools/gen_tables.py and serve the repository root ' +
      '(python3 -m http.server 8000).');
  }
  const meta = await metaRes.json();
  if (!meta.odd || !meta.radial || !meta.momentum) {
    throw new Error('orbital_tables.json lacks the odd, radial or momentum tables: it is out of date ' +
      'or cached. Hard-reload the page, or run python3 tools/gen_tables.py.');
  }
  const expected = meta.rows.length * meta.N * 4;
  const data = toFloats(await binRes.arrayBuffer(), 'orbital_tables.bin', expected);

  const [odd, radial, momentum] = await Promise.all([
    [meta.odd.file, expected], [meta.radial.file, expected / 4], [meta.momentum.file, expected / 4],
  ].map(async ([file, n]) => {
    const res = await fetch(new URL(file, base), FETCH_OPTS);
    if (!res.ok) throw new Error(`Could not load ${file} (${res.status}).`);
    return toFloats(await res.arrayBuffer(), file, n);
  }));
  return { N: meta.N, uMax: meta.uMax, rows: meta.rows, data, odd, radial, momentum };
}

// Little-endian float32 file -> Float32Array of the expected length.
function toFloats(buf, name, expected) {
  let data;
  const littleEndian = new Uint8Array(new Uint16Array([1]).buffer)[0] === 1;
  if (littleEndian) {
    data = new Float32Array(buf);
  } else {
    const view = new DataView(buf);
    data = new Float32Array(buf.byteLength / 4);
    for (let i = 0; i < data.length; i++) data[i] = view.getFloat32(4 * i, true);
  }
  if (data.length !== expected) {
    throw new Error(`${name} has ${data.length} floats, expected ${expected}`);
  }
  return data;
}
