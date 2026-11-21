// Minimal MP4 (ISO BMFF) writer for one H.264 video track.
//
// Input is what WebCodecs' VideoEncoder produces with avc: { format: 'avc' }:
// length-prefixed access units plus the avcC record
// (decoderConfig.description). The file is non-fragmented with the moov box
// before mdat ("fast start"), one sample per chunk, constant frame duration.
// If the encoder emits frames out of presentation order (B-frames), a ctts
// box carries the composition offsets.
//
// Box layout:
//   ftyp
//   moov: mvhd, trak: tkhd, mdia: mdhd, hdlr,
//         minf: vmhd, dinf/dref, stbl: stsd(avc1/avcC), stts, [ctts], stss, stsc, stsz, stco
//   mdat

const enc = new TextEncoder();

// Concatenate byte arrays.
function concat(parts) {
  const n = parts.reduce((s, p) => s + p.length, 0);
  const out = new Uint8Array(n);
  let o = 0;
  for (const p of parts) { out.set(p, o); o += p.length; }
  return out;
}

// Big-endian fields: [bits, value] pairs (8, 16, 32 bit unsigned or signed).
function fields(...pairs) {
  const n = pairs.reduce((s, [bits]) => s + bits / 8, 0);
  const b = new Uint8Array(n), v = new DataView(b.buffer);
  let o = 0;
  for (const [bits, x] of pairs) {
    if (bits === 8) v.setUint8(o, x);
    else if (bits === 16) v.setUint16(o, x);
    else v.setUint32(o, x >>> 0);
    o += bits / 8;
  }
  return b;
}

const u8 = x => [8, x], u16 = x => [16, x], u32 = x => [32, x];

function box(type, ...payload) {
  const body = concat(payload);
  return concat([fields(u32(8 + body.length)), enc.encode(type), body]);
}

// Full box: version and 24-bit flags before the payload.
function fullBox(type, version, flags, ...payload) {
  return box(type, fields(u8(version), u8(flags >> 16), u16(flags & 0xffff)), ...payload);
}

const MATRIX = fields(u32(0x10000), u32(0), u32(0), u32(0), u32(0x10000), u32(0), u32(0), u32(0), u32(0x40000000));

/**
 * Build an MP4 file.
 *   width, height   coded size
 *   timescale       ticks per second (e.g. fps * 1000)
 *   frameDuration   ticks per frame
 *   description     avcC bytes (ArrayBuffer or view) from the encoder
 *   chunks          [{ data: Uint8Array, pts: ticks, key: boolean }] in decode order
 * Returns a Uint8Array.
 */
export function muxMp4({ width, height, timescale, frameDuration, description, chunks }) {
  if (!chunks.length) throw new Error('muxMp4: no frames');
  const avcC = new Uint8Array(description instanceof ArrayBuffer ? description : description.buffer,
    description.byteOffset || 0, description.byteLength);
  const n = chunks.length;
  const duration = n * frameDuration;

  // Decode timestamps are evenly spaced; composition offsets only if needed.
  const minPts = Math.min(...chunks.map(c => c.pts));
  const offsets = chunks.map((c, i) => c.pts - minPts - i * frameDuration);
  const needCtts = offsets.some(o => o !== 0);

  const sampleEntry = box('avc1',
    fields(u32(0), u16(0), u16(1)),                     // reserved[6], data_reference_index
    fields(u16(0), u16(0), u32(0), u32(0), u32(0)),     // pre_defined, reserved, pre_defined[3]
    fields(u16(width), u16(height), u32(0x00480000), u32(0x00480000), u32(0), u16(1)),
    new Uint8Array(32),                                 // compressorname
    fields(u16(0x18), u16(0xffff)),                     // depth, pre_defined = -1
    box('avcC', avcC));

  const stts = fullBox('stts', 0, 0, fields(u32(1), u32(n), u32(frameDuration)));
  const ctts = needCtts
    ? fullBox('ctts', 1, 0, fields(u32(n)), concat(offsets.map(o => fields(u32(1), u32(o)))))
    : new Uint8Array(0);
  const keys = chunks.map((c, i) => (c.key ? i + 1 : 0)).filter(Boolean);
  const stss = fullBox('stss', 0, 0, fields(u32(keys.length)), concat(keys.map(k => fields(u32(k)))));
  const stsc = fullBox('stsc', 0, 0, fields(u32(1), u32(1), u32(1), u32(1)));
  const stsz = fullBox('stsz', 0, 0, fields(u32(0), u32(n)), concat(chunks.map(c => fields(u32(c.data.length)))));

  // stco needs absolute file offsets, which depend on the moov size; its
  // size does not depend on the values, so build moov once with zeros to
  // measure it, then again with the real offsets.
  const ftyp = box('ftyp', enc.encode('isom'), fields(u32(0x200)), enc.encode('isomiso2avc1mp41'));
  const moovWith = stcoOffsets => box('moov',
    fullBox('mvhd', 0, 0, fields(u32(0), u32(0), u32(timescale), u32(duration), u32(0x10000), u16(0x100), u16(0)),
      new Uint8Array(8), MATRIX, new Uint8Array(24), fields(u32(2))),
    box('trak',
      fullBox('tkhd', 0, 3, fields(u32(0), u32(0), u32(1), u32(0), u32(duration)),
        new Uint8Array(8), fields(u16(0), u16(0), u16(0), u16(0)), MATRIX,
        fields(u32(width * 65536), u32(height * 65536))),
      box('mdia',
        fullBox('mdhd', 0, 0, fields(u32(0), u32(0), u32(timescale), u32(duration), u16(0x55c4), u16(0))),
        fullBox('hdlr', 0, 0, fields(u32(0)), enc.encode('vide'), new Uint8Array(12), enc.encode('VideoHandler\0')),
        box('minf',
          fullBox('vmhd', 0, 1, fields(u16(0), u16(0), u16(0), u16(0))),
          box('dinf', fullBox('dref', 0, 0, fields(u32(1)), fullBox('url ', 0, 1))),
          box('stbl',
            fullBox('stsd', 0, 0, fields(u32(1)), sampleEntry),
            stts, ctts, stss, stsc, stsz,
            fullBox('stco', 0, 0, fields(u32(n)), concat(stcoOffsets.map(o => fields(u32(o))))))))));

  const moovSize = moovWith(new Array(n).fill(0)).length;
  const mdatSize = 8 + chunks.reduce((s, c) => s + c.data.length, 0);
  let o = ftyp.length + moovSize + 8;
  if (o + mdatSize > 0xffffffff) throw new Error('muxMp4: file larger than 4 GB');
  const stco = chunks.map(c => { const at = o; o += c.data.length; return at; });

  return concat([ftyp, moovWith(stco), fields(u32(mdatSize)), enc.encode('mdat'), ...chunks.map(c => c.data)]);
}
