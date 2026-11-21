// Loads the moment tables written by tools/gen_tables.py:
//   data/orbital_tables.json      metadata (N, uMax, rows, channels, odd)
//   data/orbital_tables.bin       float32 little-endian, [row][i][RGBA]
//   data/orbital_tables_odd.bin   same layout, odd half-line moments
// Texels hold (I0, I2, I4, I6) and (J1, J3, J5, 0) of one radial function
// at u = i * uMax / (N - 1).

const DEFAULT_BASE = new URL('../data/', import.meta.url);

export async function loadTables(base = DEFAULT_BASE) {
  const [metaRes, binRes] = await Promise.all([
    fetch(new URL('orbital_tables.json', base)),
    fetch(new URL('orbital_tables.bin', base)),
  ]);
  if (!metaRes.ok || !binRes.ok) {
    throw new Error(`Could not load orbital tables (${metaRes.status}/${binRes.status}). ` +
      'Run python3 tools/gen_tables.py and serve the repository root ' +
      '(python3 -m http.server 8000).');
  }
  const meta = await metaRes.json();
  const expected = meta.rows.length * meta.N * 4;
  const data = toFloats(await binRes.arrayBuffer(), 'orbital_tables.bin', expected);

  const oddRes = await fetch(new URL(meta.odd.file, base));
  if (!oddRes.ok) throw new Error(`Could not load ${meta.odd.file} (${oddRes.status}).`);
  const odd = toFloats(await oddRes.arrayBuffer(), meta.odd.file, expected);
  return { N: meta.N, uMax: meta.uMax, rows: meta.rows, data, odd };
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
