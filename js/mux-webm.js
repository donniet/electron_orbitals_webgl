// Minimal WebM (Matroska) writer for one VP8 / VP9 / AV1 video track.
//
// Input is what WebCodecs' VideoEncoder produces: one encoded frame per
// chunk (plus, for AV1, the av1C record as decoderConfig.description).
// Layout, all sizes known (the file is written in one go):
//   EBML header
//   Segment: SeekHead (-> Info, Tracks, Cues), Info, Tracks,
//            Cluster per keyframe (Timestamp, SimpleBlocks), Cues
// Timestamps are in milliseconds (TimestampScale 1e6 ns), the WebM norm.

const enc = new TextEncoder();

function concat(parts) {
  const n = parts.reduce((s, p) => s + p.length, 0);
  const out = new Uint8Array(n);
  let o = 0;
  for (const p of parts) { out.set(p, o); o += p.length; }
  return out;
}

// Element size as an EBML variable-length integer (shortest form).
function vint(n) {
  for (let len = 1; len <= 8; len++) {
    if (n < 2 ** (7 * len) - 1) {
      const b = new Uint8Array(len);
      let x = n;
      for (let i = len - 1; i >= 0; i--) { b[i] = x % 256; x = Math.floor(x / 256); }
      b[0] |= 1 << (8 - len);
      return b;
    }
  }
  throw new Error('EBML size too large');
}

// Unsigned integer, big-endian, in `width` bytes (default: shortest).
function uint(x, width = 0) {
  const bytes = [];
  do { bytes.unshift(x % 256); x = Math.floor(x / 256); } while (x > 0);
  while (bytes.length < width) bytes.unshift(0);
  return new Uint8Array(bytes);
}

function float64(x) {
  const b = new Uint8Array(8);
  new DataView(b.buffer).setFloat64(0, x);
  return b;
}

// Element ids are written with their length marker included, as in the spec.
function el(id, ...payload) {
  const body = concat(payload.map(p => (typeof p === 'string' ? enc.encode(p) : p)));
  return concat([uint(id), vint(body.length), body]);
}

const ID = {
  EBML: 0x1A45DFA3, EBMLVersion: 0x4286, EBMLReadVersion: 0x42F7, EBMLMaxIDLength: 0x42F2,
  EBMLMaxSizeLength: 0x42F3, DocType: 0x4282, DocTypeVersion: 0x4287, DocTypeReadVersion: 0x4285,
  Segment: 0x18538067, SeekHead: 0x114D9B74, Seek: 0x4DBB, SeekID: 0x53AB, SeekPosition: 0x53AC,
  Info: 0x1549A966, TimestampScale: 0x2AD7B1, Duration: 0x4489, MuxingApp: 0x4D80, WritingApp: 0x5741,
  Tracks: 0x1654AE6B, TrackEntry: 0xAE, TrackNumber: 0xD7, TrackUID: 0x73C5, TrackType: 0x83,
  FlagLacing: 0x9C, DefaultDuration: 0x23E383, CodecID: 0x86, CodecPrivate: 0x63A2,
  Video: 0xE0, PixelWidth: 0xB0, PixelHeight: 0xBA,
  Cluster: 0x1F43B675, Timestamp: 0xE7, SimpleBlock: 0xA3,
  Cues: 0x1C53BB6B, CuePoint: 0xBB, CueTime: 0xB3, CueTrackPositions: 0xB7, CueTrack: 0xF7,
  CueClusterPosition: 0xF1,
};

export const WEBM_CODEC_IDS = { vp8: 'V_VP8', vp09: 'V_VP9', av01: 'V_AV1' };

/**
 * Build a WebM file.
 *   width, height  frame size
 *   fps            frames per second (for durations)
 *   codec          WebCodecs codec string ('vp8', 'vp09...', 'av01...')
 *   description    optional codec private data (AV1: av1C)
 *   chunks         [{ data: Uint8Array, pts: microseconds, key: boolean }] in order
 * Returns a Uint8Array.
 */
export function muxWebm({ width, height, fps, codec, description, chunks }) {
  if (!chunks.length) throw new Error('muxWebm: no frames');
  const codecId = WEBM_CODEC_IDS[codec.split('.')[0]];
  if (!codecId) throw new Error(`muxWebm: unsupported codec ${codec}`);
  const ms = us => Math.round(us / 1000);
  const durationMs = ms(chunks[chunks.length - 1].pts) + 1000 / fps;

  const header = el(ID.EBML,
    el(ID.EBMLVersion, uint(1)), el(ID.EBMLReadVersion, uint(1)),
    el(ID.EBMLMaxIDLength, uint(4)), el(ID.EBMLMaxSizeLength, uint(8)),
    el(ID.DocType, 'webm'), el(ID.DocTypeVersion, uint(4)), el(ID.DocTypeReadVersion, uint(2)));

  const info = el(ID.Info,
    el(ID.TimestampScale, uint(1000000)), el(ID.Duration, float64(durationMs)),
    el(ID.MuxingApp, 'electron-orbitals'), el(ID.WritingApp, 'electron-orbitals'));

  const priv = description
    ? [el(ID.CodecPrivate, new Uint8Array(description instanceof ArrayBuffer ? description : description.buffer,
        description.byteOffset || 0, description.byteLength))]
    : [];
  const tracks = el(ID.Tracks, el(ID.TrackEntry,
    el(ID.TrackNumber, uint(1)), el(ID.TrackUID, uint(1)), el(ID.TrackType, uint(1)),
    el(ID.FlagLacing, uint(0)), el(ID.DefaultDuration, uint(Math.round(1e9 / fps))),
    el(ID.CodecID, codecId), ...priv,
    el(ID.Video, el(ID.PixelWidth, uint(width)), el(ID.PixelHeight, uint(height)))));

  // Clusters: a new one at every keyframe (and before the 16-bit relative
  // block timestamp would overflow).
  const clusters = [];
  let cur = null;
  for (const c of chunks) {
    const t = ms(c.pts);
    if (!cur || c.key || t - cur.t > 30000) { cur = { t, key: c.key, blocks: [] }; clusters.push(cur); }
    const rel = t - cur.t;
    cur.blocks.push(el(ID.SimpleBlock, vint(1), new Uint8Array([(rel >> 8) & 0xff, rel & 0xff, c.key ? 0x80 : 0]), c.data));
  }
  const clusterBytes = clusters.map(cl => el(ID.Cluster, el(ID.Timestamp, uint(cl.t)), ...cl.blocks));

  // Positions are relative to the start of the Segment's payload. The
  // SeekHead uses fixed 8-byte positions so its size is known up front.
  const seekHeadWith = (pInfo, pTracks, pCues) => el(ID.SeekHead,
    ...[[ID.Info, pInfo], [ID.Tracks, pTracks], [ID.Cues, pCues]].map(([id, pos]) =>
      el(ID.Seek, el(ID.SeekID, uint(id)), el(ID.SeekPosition, uint(pos, 8)))));
  const seekSize = seekHeadWith(0, 0, 0).length;
  const pInfo = seekSize, pTracks = pInfo + info.length;
  let pos = pTracks + tracks.length;
  const cuePoints = [];
  clusters.forEach((cl, i) => {
    if (cl.key) {
      cuePoints.push(el(ID.CuePoint, el(ID.CueTime, uint(cl.t)),
        el(ID.CueTrackPositions, el(ID.CueTrack, uint(1)), el(ID.CueClusterPosition, uint(pos)))));
    }
    pos += clusterBytes[i].length;
  });
  const cues = el(ID.Cues, ...cuePoints);
  const segmentBody = [seekHeadWith(pInfo, pTracks, pos), info, tracks, ...clusterBytes, cues];
  const segmentSize = segmentBody.reduce((s, p) => s + p.length, 0);
  return concat([header, uint(ID.Segment), vint(segmentSize), ...segmentBody]);
}
