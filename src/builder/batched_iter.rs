/// Iterate `packed_bytes` over factors of `k`, yielding 8 consecutive k-mers
/// (one batch) at a time as a plain `[u32; 8]`.
pub struct BatchedKmerIterator<'a> {
    packed_bytes: &'a [u8],
    pos_bytes: usize,
    bit_offset: u32, // Tracks the constant intra-byte shift
    filter_mask: u32,
}

impl<'a> BatchedKmerIterator<'a> {
    pub const fn new(k: usize, packed_bytes: &'a [u8], initial_shift_bases: usize) -> Self {
        let filter_mask = u32::MAX >> (32 - 2 * k);
        let pos_bytes = initial_shift_bases / 4;
        let bit_offset = ((initial_shift_bases % 4) * 2) as u32;

        Self {
            packed_bytes,
            pos_bytes,
            bit_offset,
            filter_mask,
        }
    }
}

impl<'a> Iterator for BatchedKmerIterator<'a> {
    type Item = [u32; 8];

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos_bytes + 8 > self.packed_bytes.len() {
            return None;
        }

        let mut window = unsafe {
            self.packed_bytes
                .as_ptr()
                .add(self.pos_bytes)
                .cast::<u64>()
                .read_unaligned()
        };
        window >>= self.bit_offset;

        let m = self.filter_mask;
        // Extract the 8 kmers
        let kmers = [
            window as u32 & m,
            (window >> 2) as u32 & m,
            (window >> 4) as u32 & m,
            (window >> 6) as u32 & m,
            (window >> 8) as u32 & m,
            (window >> 10) as u32 & m,
            (window >> 12) as u32 & m,
            (window >> 14) as u32 & m,
        ];

        self.pos_bytes += 2;

        Some(kmers)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::utils::kmer_to_ascii;
    use bytemuck::cast_slice;
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
            let (packed_data, _) = parser.get_dna_packed().bits();
            let packed_bytes = cast_slice::<u128, u8>(packed_data);
            let mut kmer_it = BatchedKmerIterator::new(k, packed_bytes, 0);

            let all_kmers: Vec<u32> = kmer_it.by_ref().flatten().collect();

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
    fn batched_prefix_iterator() {
        let fasta: &[u8] = b">test_seq\nCCCTGAGTACGGAAAGCGCGAACGCAGATGCCCTATCGATACGTGGCAAGAGTGTTGTCCAAAGGGGCTACGCCCCTATTGAGTATTTACTATTGATTGTTAGATGTGAGTGCGTCTCAATCCTGCCGTTACTTGACCGTTTATGAGTTTATTATAGTCGTTAATATCTGGTCGAGACGGTGTAAAATACGCTATGCGACACCTGTCGTATATCAGAGAAAAGGGTCGATTCTCAATAAATATCGCCCTCTAAACCAGTTTAGGATGCTCTGGAGCCGAAGGATGGGTTCTTGCAGAATACATCACTTCTAGTAAGCGTCAGGCAAACGGCTTTAACCACCTTAGAAAGGGGCAATCACCCAAAGAATACAGTTGAGTAACGATTGTAAAAATAATGTAACAATGCATCAGTAGGAATCACCTTCACTTTCTTTGTATAGGAGTACGCACTCTTGTGGATACCCTCCGAACTACATACACGGTCCCAGTAACAGAGCT";
        for k in 1..=16 {
            check_kmer_iterator(fasta, k);
        }
    }
}
