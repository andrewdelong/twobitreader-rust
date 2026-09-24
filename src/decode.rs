// Standard library
use std::iter::zip;
use std::mem::MaybeUninit;

// Dependencies
use seq_macro::seq;

// Constants related to decoding.
pub(crate) const BITS_PER_U8: usize = 8;
pub(crate) const BITS_PER_U32: usize = 32;
pub(crate) const BITS_PER_NUC: usize = 2;
pub(crate) const NUCS_PER_U8: usize = BITS_PER_U8 / BITS_PER_NUC;

// Lookup table to decode two bits into an ASCII character: T=00 C=01 A=10 G=11.
const DECODE_U2: [u8; 4] = *b"TCAG";

// Lookup table to decode eight bits into four ASCII characters at once.
// For example, DECODE_U8[0b_11_00_00_10] gives a u32 representing 'GTTA'.
// The seq! macro expands at compile time to:
//   const DECODE_U8: [u32; 256] = [                     // as binary:
//       u32::from_le_bytes([b'T', b'T', b'T', b'T']),   // 00 00 00 00
//       u32::from_le_bytes([b'T', b'T', b'T', b'C']),   // 00 00 00 01
//       u32::from_le_bytes([b'T', b'T', b'T', b'A']),   // 00 00 00 10
//       u32::from_le_bytes([b'T', b'T', b'T', b'G']),   // 00 00 00 11
//       u32::from_le_bytes([b'T', b'T', b'C', b'T']),   // 00 00 01 00
//       ...
//       u32::from_le_bytes([b'G', b'G', b'A', b'G']),   // 11 11 10 11
//       u32::from_le_bytes([b'G', b'G', b'G', b'T']),   // 11 11 11 00
//       u32::from_le_bytes([b'G', b'G', b'G', b'C']),   // 11 11 11 01
//       u32::from_le_bytes([b'G', b'G', b'G', b'A']),   // 11 11 11 10
//       u32::from_le_bytes([b'G', b'G', b'G', b'G']),   // 11 11 11 11
//   ];
#[allow(clippy::identity_op)] // Omit clippy false positives (confused by seq! macro)
#[allow(clippy::erasing_op)]
#[allow(clippy::eq_op)]
const DECODE_U8: [u32; 256] = seq!(i in 0..256 {[#(
    u32::from_le_bytes([DECODE_U2[(i / 64) % 4],
                        DECODE_U2[(i / 16) % 4],
                        DECODE_U2[(i / 4)  % 4],
                        DECODE_U2[(i / 1)  % 4]]),
)*]});

// Decode 2-bit packed dna, starting at the given position.
// The contents of dst are replaced by ASCII nucleotides.
pub(crate) fn decode(start: usize, dna: &[u8], dst: &mut [MaybeUninit<u8>]) {
    // Split dst into quads (groups of four nucleotides), plus any trailing u8s left over.
    // Quads are the unit the main decoding loop works in, because one packed source byte holds
    // exactly four nucleotides, so each iteration reads one byte and writes four.
    let (dst_quads, dst_tail) = dst.as_chunks_mut::<NUCS_PER_U8>();

    // Calculate the start byte / start alignment corresponding to the output slices.
    // - start_byte: index of first byte to be decoded
    // - start_align: index of first nucleotide within byte
    // - tail_start_byte: index of first byte to be decoded after quads
    //
    // Every quad begins NUCS_PER_U8 nucleotides after the one before it, so they all share the
    // same alignment within their source byte, namely start_align. The tail does too.
    let start_byte = start / NUCS_PER_U8;
    let start_align = start % NUCS_PER_U8;
    let tail_start_byte = start_byte + dst_quads.len();

    // Decode as many quads as possible.
    // - The decoding loop is specialized for each of the four possible start_align values {0,1,2,3},
    //   which allows the compiler to bit-shift by a known constant, rather than by a variable amount.
    match start_align {
        0 => decode_quads::<0>(start_byte, dna, dst_quads),
        1 => decode_quads::<1>(start_byte, dna, dst_quads),
        2 => decode_quads::<2>(start_byte, dna, dst_quads),
        3 => decode_quads::<3>(start_byte, dna, dst_quads),
        _ => unreachable!(),
    };

    // Decode whatever trailing nucleotides remained to be filled after quads.
    decode_subquad(tail_start_byte, start_align, dna, dst_tail);
}

// Decodes up to 3 nucleotides (i.e., less than a quad) from 2-bit packed DNA.
// - start_byte: index of byte within dna to decode.
// - start_align: index of start nucleotide within the start byte.
// - dna: packed bytes of DNA for a 2bit sequence record.
// - dst: destination to decode ASCII nucleotides into.
fn decode_subquad(start_byte: usize, start_align: usize, dna: &[u8], dst: &mut [MaybeUninit<u8>]) {
    if dst.is_empty() {
        return;
    }

    // Indicate that final loop can be unrolled; assume! macro would work fine too.
    assert!(dst.len() < NUCS_PER_U8);

    // Decode the start byte as a quad
    let mut nucs = DECODE_U8[dna[start_byte] as usize];

    if start_align > 0 {
        // Shift out any unwanted leading nucleotides from the low-order bits of nucs.
        nucs >>= BITS_PER_U8 * start_align;

        // If not enough nucleotides remain to fill dst, incorporate more from the next byte.
        if dst.len() > NUCS_PER_U8 - start_align {
            let extra_nucs = DECODE_U8[dna[start_byte + 1] as usize];
            nucs |= extra_nucs << (BITS_PER_U32 - BITS_PER_U8 * start_align);
        }
    }

    // Write the relevant ASCII bytes from nucs into dst.
    for d in dst {
        d.write((nucs & 0xff) as u8);
        nucs >>= BITS_PER_U8;
    }
}

// Decode a span of 2-bit encoded nucleotide sequence that spans multiple bytes in memory.
// - start_byte: index of byte within dna to decode.
// - START_ALIGN: index of start nucleotide within the start byte.
// - dna: packed bytes of DNA for a 2bit sequence record.
// - dst: destination to decode quads of ASCII nucleotides into.
fn decode_quads<const START_ALIGN: usize>(start_byte: usize, dna: &[u8], dst: &mut [[MaybeUninit<u8>; NUCS_PER_U8]]) {
    if dst.is_empty() {
        return;
    }

    // Slice all bytes of packed dna that will be needed (extra byte if start not byte-aligned)
    let src_len = dst.len() + START_ALIGN.min(1);
    let src = &dna[start_byte..start_byte + src_len];

    if START_ALIGN == 0 {
        // Source quads from dna are all byte-aligned, so just decode directly into dst.
        for (s, d) in zip(src, dst) {
            d.write_copy_of_slice(&DECODE_U8[*s as usize].to_le_bytes());
        }
    } else {
        // Source quads from dna are not byte-aligned, so each decode straddles two src bytes.
        // Each output quad is produced by combining the bits of the two most recently decoded quads.
        let mut prev = DECODE_U8[src[0] as usize];
        for (s, d) in zip(&src[1..], dst) {
            let next = DECODE_U8[*s as usize];
            d.write_copy_of_slice(&combine_quads::<START_ALIGN>(prev, next).to_le_bytes());
            prev = next;
        }
    }
}

// Given two quads with prev bytes 0..3 and next bytes 4..7 this function will
// extract four consecutive bytes starting at START_ALIGN. For example,
// if START_ALIGN=2 then the returned quad is formed from bytes 2..5.
#[inline]
fn combine_quads<const START_ALIGN: usize>(prev: u32, next: u32) -> u32 {
    // Emits a single shrd instruction on optimized x86 builds.
    (prev >> (BITS_PER_U8 * START_ALIGN)) | (next << (BITS_PER_U32 - BITS_PER_U8 * START_ALIGN))
}
