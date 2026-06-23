use bytemuck::cast_slice;
use wide::u32x8;

/// Iterate `packed_bytes` over factors of k `prefix_size`
pub struct SimdKmerIterator<'a> {
    packed_bytes: &'a [u8],
    pos_bytes: usize,
    bit_offset: u32, // Tracks the constant intra-byte shift
    filter_mask: u32x8,
}

impl<'a> SimdKmerIterator<'a> {
    pub fn new(k: usize, packed_data: &'a [u128], initial_shift_bases: usize) -> Self {
        let packed_bytes = cast_slice::<u128, u8>(packed_data);
        let filter_prefix_mask = u32::MAX >> (32 - 2 * k);
        let filter_mask = u32x8::from([filter_prefix_mask; 8]);
        let pos_bytes = initial_shift_bases / 4;
        let bit_offset = ((initial_shift_bases % 4) * 2) as u32;

        return Self {
            packed_bytes,
            pos_bytes,
            bit_offset,
            filter_mask,
        }
    }

    pub fn remainder_bytes(&self) -> &'a [u8] {
        &self.packed_bytes[self.pos_bytes..]
    }
}

impl<'a> Iterator for SimdKmerIterator<'a> {
    type Item = u32x8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos_bytes + 8 > self.packed_bytes.len() {
            return None;
        }

        // Create current u64 window
        let window_bytes: [u8; 8] = self.packed_bytes[self.pos_bytes..self.pos_bytes + 8]
            .try_into()
            .unwrap();

        let mut window = u64::from_le_bytes(window_bytes);
        window >>= self.bit_offset;

        // Extrtact the 8 kmers
        let kmers = [
            (window >> 0) as u32,
            (window >> 2) as u32,
            (window >> 4) as u32,
            (window >> 6) as u32,
            (window >> 8) as u32,
            (window >> 10) as u32,
            (window >> 12) as u32,
            (window >> 14) as u32,
        ];

        let v_kmers = u32x8::from(kmers);
        let result = v_kmers & self.filter_mask;

        self.pos_bytes += 2;

        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::utils::kmer_to_ascii;
    use helicase::{Config, FastxParser, HelicaseParser, ParserOptions, input::FromSlice};

    const CONFIG: Config = ParserOptions::default()
        .ignore_headers()
        .dna_string()
        .and_dna_packed()
        .config();

    fn check_kmer_iterator(fasta: &[u8], k: usize) {
        let mut parser =
            FastxParser::<CONFIG>::from_slice(fasta).expect("Failed to initialize parser");
        let mut buf = vec![0; k];

        while let Some(_) = parser.next() {
            let ascii_bytes = parser.get_dna_string();
            let (packed_bytes, _) = parser.get_dna_packed().bits();
            let mut kmer_it = SimdKmerIterator::new(k, packed_bytes, 0);

            let all_kmers: Vec<u32> = kmer_it
                        .by_ref()
                        .flat_map(|v_kmer| v_kmer.to_array())
                        .collect();

            ascii_bytes
                .windows(k)
                .zip(all_kmers)
                .for_each(|(ascii, kmer)| {
                    kmer_to_ascii(kmer, k, &mut buf);
                    assert_eq!(ascii, &buf[..k])
                });
        }
    }

    #[test]
    fn simd_prefix_iterator() {
        let fasta: &[u8] = b">test_seq\nCCCTGAGTACGGAAAGCGCGAACGCAGATGCCCTATCGATACGTGGCAAGAGTGTTGTCCAAAGGGGCTACGCCCCTATTGAGTATTTACTATTGATTGTTAGATGTGAGTGCGTCTCAATCCTGCCGTTACTTGACCGTTTATGAGTTTATTATAGTCGTTAATATCTGGTCGAGACGGTGTAAAATACGCTATGCGACACCTGTCGTATATCAGAGAAAAGGGTCGATTCTCAATAAATATCGCCCTCTAAACCAGTTTAGGATGCTCTGGAGCCGAAGGATGGGTTCTTGCAGAATACATCACTTCTAGTAAGCGTCAGGCAAACGGCTTTAACCACCTTAGAAAGGGGCAATCACCCAAAGAATACAGTTGAGTAACGATTGTAAAAATAATGTAACAATGCATCAGTAGGAATCACCTTCACTTTCTTTGTATAGGAGTACGCACTCTTGTGGATACCCTCCGAACTACATACACGGTCCCAGTAACAGAGCT";
        for k in 1..=16 {
            check_kmer_iterator(fasta, k);
        }
    }
}
