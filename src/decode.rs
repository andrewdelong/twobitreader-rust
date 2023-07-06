// Standard library
use std::iter::zip;

// Dependencies
use seq_macro::seq;

// Constants related to decoding.
pub const BITS_PER_U8: usize = 8;
pub const BITS_PER_U32: usize = 32;
pub const BITS_PER_NUC: usize = 2;
pub const NUCS_PER_U8: usize = BITS_PER_U8 / BITS_PER_NUC;

// Lookup table to decode two bits into an ASCII character: T=00 C=01 A=10 G=11.
const DECODE_U2: [u8; 4] = [b'T', b'C', b'A', b'G'];

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
const DECODE_U8: [u32; 256] = seq!(i in 0..256 {[#(
    u32::from_le_bytes([DECODE_U2[i / 64 % 4],
                        DECODE_U2[i / 16 % 4],
                        DECODE_U2[i / 4  % 4],
                        DECODE_U2[i / 1  % 4]]),
)*]});

// Decode 2-bit packed dna, starting at the given position.
// The contents of dst are replaced by ASCII nucleotides.
pub fn decode(start: usize, dna: &[u8], dst: &mut[u8]) {

    // Transmute dst into a u32 slice (quads), plus any leading/trailing u8s needed to align it.
    // This allows the main quad-based decoding loop to perform correctly-aligned writes,
    // which is important for speed and for portability.
    let (dst_prequads, dst_quads, dst_postquads) = unsafe { dst.align_to_mut() };
    
    // Calculate the start byte / start alignment corresponding each of the output slices.
    let start_byte = start / NUCS_PER_U8;               // Index of first byte to be decoded
    let start_align = start % NUCS_PER_U8;              // Index of first nucleotide within byte
    let quads_start = start + dst_prequads.len();       // Position of first nucleotide of u32-aligned dst
    let quads_start_byte = quads_start / NUCS_PER_U8;   // Index of first byte to be decoded as quad
    let quads_start_align = quads_start % NUCS_PER_U8;  // Index of first nucleotide within byte to be decoded as quad
    let postquads_start_byte = quads_start_byte + dst_quads.len();  // Index of first byte to be decoded after quads

    // Decode enough just enough leading nucleotides to bring dst into u32 alignment.
    decode_subquad(start_byte, start_align, dna, dst_prequads);

    // Decode as many quads as possible.
    // - The decoding loop is specialized for each of the four possible start_align values {0,1,2,3},
    //   which allows the compiler to bit-shift by a known constant, rather than by a variable amount.
    match quads_start_align {
        0 => decode_quads::<0>(quads_start_byte, dna, dst_quads),
        1 => decode_quads::<1>(quads_start_byte, dna, dst_quads),
        2 => decode_quads::<2>(quads_start_byte, dna, dst_quads),
        3 => decode_quads::<3>(quads_start_byte, dna, dst_quads),
        _ => unreachable!(),
    };

    // Decode whatever trailing nucleotides remained to be filled after quads.
    decode_subquad(postquads_start_byte, quads_start_align, dna, dst_postquads);
}

// Decodes up to 3 nucleotides (i.e., less than a quad) from 2-bit packed DNA.
// - start_byte: index of byte within dna to decode.
// - start_align: index of start nucleotide within the start byte.
// - dna: packed bytes of DNA for a 2bit sequence record.
// - dst: destination to decode ASCII nucleotides into.
fn decode_subquad(start_byte: usize, start_align: usize, dna: &[u8], dst: &mut [u8]) {
    if dst.len() == 0 { return; }
    assert!(dst.len() < NUCS_PER_U8);

    // Decode the start byte as a quad
    let mut nucs = DECODE_U8[dna[start_byte] as usize];

    if start_align > 0 {
        // Shift out any unwanted leading nucleotides from the low-order bits of nucs.
        nucs >>= BITS_PER_U8 * start_align;

        // If not enough nucleotides remain to fill dst, incorporate more from the next byte.
        if dst.len() > NUCS_PER_U8 - start_align {
            let extra_nucs = DECODE_U8[dna[start_byte+1] as usize];
            nucs |= extra_nucs << (BITS_PER_U32 - BITS_PER_U8 * start_align);
        }
    }

    // Write the relevant ASCII bytes from nucs into dst.
    for d in dst {
        *d =  (nucs & 0xff) as u8;
        nucs >>= BITS_PER_U8;
    }
}

// Decode a span of 2-bit encoded nucleotide sequence that spans multiple bytes in memory.
// - start_byte: index of byte within dna to decode.
// - START_ALIGN: index of start nucleotide within the start byte.
// - dna: packed bytes of DNA for a 2bit sequence record.
// - dst: destination to decode quads of ASCII nucleotides into.
fn decode_quads<const START_ALIGN: usize>(start_byte: usize, dna: &[u8], dst: &mut [u32]) {
    if dst.len() == 0 { return; }

    // Slice all bytes of packed dna that will be needed (extra byte if start not byte-aligned)
    let src_len = dst.len() + START_ALIGN.min(1);
    let src = &dna[start_byte..start_byte+src_len];

    if START_ALIGN == 0 {
        // Source quads from dna are all byte-aligned, so just decode directly into dst.
        for (s, d) in zip(src, dst) {
            *d = DECODE_U8[*s as usize];
        }
    } else {
        // Source quads from dna are not byte-aligned, so each decode straddles two src bytes.
        // Each output quad is produced by combining the bits of the two most recently decoded quads.
        let mut prev = DECODE_U8[src[0] as usize];
        for (s, d) in zip(&src[1..], dst) {
            let next = DECODE_U8[*s as usize];
            *d = combine_quads::<START_ALIGN>(prev, next);
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
    (prev >> BITS_PER_U8 * START_ALIGN) | (next << (BITS_PER_U32 - BITS_PER_U8 * START_ALIGN))
}



        // Example 1: aligned start and aligned end on little-endian architecture
        //
        //   src: ATTG ATCC                  src4.len = 2
        //        0000 1111              <-- index of src byte
        //        0123 4567              <-- index of 2-bit code
        //
        //   dst: ???? ????                  dst4.len = 2
        //        0123 4567              <-- index of dst byte
        //
        // Step 1:
        //   x = DECODE[src[0]]
        //   dst[0:4] = x
        //   dst: ATTG ????
        //        0123 4567
        //
        // Step 2:
        //   x = DECODE[src[1]]
        //   dst[4:8] = x
        //   dst: ATTG ATCC
        //        0123 4567
        //
        //
        // Example 2: aligned start and unaligned end on little-endian architecture
        //
        //   src: ATTG ATCC GCCA              src4.len = 2
        //        0000 1111 2222          <-- index of src byte
        //        0123 4567 89A           <-- index of 2-bit code     divrem(b,4) = (2,3)  divrem(b+3,4) = (3,2)   "decode-and-shift 2 times, then 3 extra nucleotides"
        //
        //   dst: ???? ???? ???               dst4.len = 2
        //        0123 4567 89A           <-- index of dst byte
        //
        // Step 1:
        //   x = DECODE[src[0]]           <-- ATTG
        //   dst[0:4] = x
        //   dst: ATTG ???? ???
        //        0123 4567 89A
        //
        // Step 2:
        //   x = DECODE[src[1]]           <-- ATCC
        //   dst[4:8] = x
        //   dst: ATTG ATCC ???
        //        0123 4567 89A
        //
        // Step 3:
        //   x = DECODE[src[2]]           <-- GCCA
        //   dst[8] = x >> 0 & 0xff       <-- G
        //   dst[9] = x >> 8 & 0xff       <-- C
        //   dst[10] = x >> 16 & 0xff     <-- C
        //   dst: ATTG ATCC GCC
        //        0123 4567 89A
        //
        //
        // Example 3: unaligned start and aligned end on little-endian architecture
        //
        //   src: ATTG ATCC GCCA          src4.len = 3   divrem(a,4) = (0,3)  divrem(b,4)  divrem(l,4) = (2,1)  => (3,0) implies need to decode 3 src bytes, since zero modulus
        //        0000 1111 2222          <-- index of src byte
        //           0 1234 5678          <-- index of 2-bit code        divrem(b,4) = (3,0)   divrem(b+3,4) = (3,3)    "decode-and-shift 3 times, then 0 extra nucleotides"
        //
        //   dst: ???? ???? ?             dst4.len = 2
        //        0123 4567 8             <-- index of dst byte
        //
        // Step 0:
        //   p = 24                       <-- amount to shift 
        //   x = DECODE[src[0]]           <-- ATTG
        //
        // Step 1:
        //   y = DECODE[src[1]]           <-- ATCC
        //   z = shld<p>(x, y)            <-- GATC
        //   dst[0:4] = z
        //   x = y
        //   dst: GATC ???? ?
        //        0123 4567 8
        //
        // Step 2:
        //   y = DECODE[src[2]]           <-- GCCA
        //   z = shld<p>(x, y)            <-- CGCC
        //   dst[4:8] = z
        //   x = y
        //   dst: GATC CGCC ?
        //        0123 4567 8
        //
        // Step 3:
        //   y = DECODE[src[3]]           <-- CTAG
        //   z = shld<p>(x, y)            <-- ACTA
        //   dst[8] = z
        //
        //
        // Example 4: unaligned start and unaligned end on little-endian architecture
        //
        //   src: ATTG ATCC GCCA CTAG     src4.len = 3   divrem(a,4) = (0,3) divrem(b,4) = (3,2) divrem(l,4) = (2,3)   => (3,2) implies need 4 src bytes, since positive modulus
        //        0000 1111 2222 3333     <-- index of src byte
        //           0 1234 5678 9A       <-- index of 2-bit code        divrem(b,4) = (3,2)   divrem(b+3,4) = (4,1)  "decode-and-shift 3 times, then two extra nucleotides"
        //
        //   dst: ???? ???? ???           dst4.len = 2
        //        0123 4567 89A           <-- index of dst byte
        //
        // Step 0:
        //   p = 24                       <-- amount to shift 
        //   x = DECODE[src[0]]           <-- ATTG
        //
        // Step 1:
        //   y = DECODE[src[1]]           <-- ATCC
        //   z = shld<p>(x, y)            <-- GATC
        //   dst[0:4] = z
        //   x = y
        //   dst: GATC ???? ???
        //        0123 4567 89A
        //
        // Step 2:
        //   y = DECODE[src[2]]           <-- GCCA
        //   z = shld<p>(x, y)            <-- CGCC
        //   dst[4:8] = z
        //   x = y
        //   dst: GATC CGCC ???
        //        0123 4567 89A
        //
        // Step 3:
        //   y = DECODE[src[3]]           <-- CTAG
        //   z = shld<p>(x, y)            <-- ACTA
        //   dst[8] = z >> 0 & 0xff       <-- A
        //   dst[9] = z >> 8 & 0xff       <-- C
        //   dst[10] = z >> 16 & 0xff     <-- T
        //
        //
        // Example 5: aligned start and unaligned within same byte
        //
        //   src: ATTG ATCC               src4.len = 0
        //        0000 1111               <-- index of src byte
        //              01                <-- index of 2-bit code
        //
        //   dst: ??                      dst4.len = 0
        //        01                      <-- index of dst byte
        //
        // Step 1:
        //   x = DECODE[src[0]]           <-- ATCC
        //   z = x << p                   <-- TCC0
        //   dst[0] = z >> 0 & 0xff       <-- T
        //   dst[1] = z >> 8 & 0xff       <-- C

        // 0000 1111 2222 3333        div(a,4) div(b,4) div(l,4) src_len src4_len (b/4-a/4)+(b%4)? 4-(a%4)  dst writes: 2 u32
        // 0123 4567                  (0,0)    (2,0)    (2,0)    2=0+2+0 2        2                0 [0123 ] [4567 ]
        //  012 3456 7                (0,1)    (2,1)    (2,0)    3=1+2+0 3        3                3 [012_3] [456_7]
        //   01 2345 67               (0,2)    (2,2)    (2,0)    3=1+2+0 3        3                2 [01_23] [45_67]
        //    0 1234 567              (0,3)    (2,3)    (2,0)    3=1+2+0 3        3                1 [0_123] [4_567]

        // 0000 1111 2222 3333        div(a,4) div(b,4) div(l,4) src_len src4_len (b/4-a/4)+(b%4)?  dst writes: 2 u32 + 1 u8
        // 0123 4567 8                (0,0)    (2,1)    (2,1)    3=0+2+1 2        3               0   [0123_] [4567_] [8]
        //  012 3456 78               (0,1)    (2,2)    (2,1)    3=1+2+0 3        3               3  [012_3] [456_7] [8]
        //   01 2345 678              (0,2)    (2,3)    (2,1)    3=1+2+0 3        3               2  [01_23] [45_67] [8]
        //    0 1234 5678             (0,3)    (3,0)    (2,1)    3=1+2+0 3        3               1  [0_123] [4_567] [8]

        // 0000 1111 2222 3333        div(a,4) div(b,4) div(l,4) src_len src4_len (b/4-a/4)+(b%4)?  dst writes: 2 u32 + 2 u8
        // 0123 4567 89               (0,0)    (2,2)    (2,2)    3=0+2+1 2        3                0 [0123_] [4567_] [8] [9 ]
        //  012 3456 789              (0,1)    (2,3)    (2,2)    3=1+2+0 3        3                3 [012_3] [456_7] [8] [9 ]
        //   01 2345 6789             (0,2)    (3,0)    (2,2)    3=1+2+0 3        3                2 [01_23] [45_67] [8] [9 ]
        //    0 1234 5678 9           (0,3)    (3,1)    (2,2)    4=1+2+1 3        4                1 [0_123] [4_567] [8] [_9]

        // 0000 1111 2222 3333        div(a,4) div(b,4) div(l,4) src_len src4_len (b/4-a/4)+(b%4)?  dst writes: 2 u32 + 3 u8
        // 0123 4567 89A              (0,0)    (2,3)    (2,3)    3=0+2+1 2        3               0  [0123_] [4567_] [8] [9] [A ]
        //  012 3456 789A             (0,1)    (3,0)    (2,3)    3=1+2+0 3        3               3  [012_3] [456_7] [8] [9] [A ]
        //   01 2345 6789 A           (0,2)    (3,1)    (2,3)    4=1+2+1 3        4               2  [01_23] [45_67] [8] [9] [_A]
        //    0 1234 5678 9A          (0,3)    (3,2)    (2,3)    4=1+2+1 3        4               1  [0_123] [4_567] [8] [_9] [A]

        // 0000 1111 2222 3333        div(a,4) div(b,4) div(l,4) src_len src4_len (b/4-a/4)+(b%4)?  dst writes: 3 u32
        // 0123 4567 89AB             (0,0)    (3,0)    (3,0)    3=0+3+0 3        3               0  [0123_] [4567_] [89AB ]
        //  012 3456 789A B           (0,1)    (3,1)    (3,0)    4=1+3+0 4        4               3  [012_3] [456_7] [89A_B]
        //   01 2345 6789 AB          (0,2)    (3,2)    (3,0)    4*1+3+0 4        4               2  [01_23] [45_67] [89_AB]
        //    0 1234 5678 9AB         (0,3)    (3,3)    (3,0)    4=1+3+0 4        4               1  [0_123] [4_567] [8_9AB]

        // Conclusions:
        //  - let dst4_len = (b-a)/4
        //  - let dst4_rem = (b-a)%4
        //  - decode4 main loop should:
        //    - read exactly dst4_len u8 values, plus one initial read if a%4 > 0
        //    - write exactly dst4_len u32 values
        //    - return u32 containing unused nucleotides
        //  - decode4 postscript should:
        //    - read 1 additional byte if dst4_rem > remaining nucs from last byte = (4-a%4)
        //    - write exactly dst4_rem u8 values
        //  - src_len is any of {dst4_len, dst4_len+1, dst4_len+2}; therefore:
        //    -  a%4? => +1 extra leading byte must be read
        //    -  b%4? => +1 extra leading byte must be read
        //    decode4 postscript should 
