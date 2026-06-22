use crate::utils::get_u64_unaligned;
use bytemuck::cast_slice;

/// Iterate `packed_bytes` over factors of k `prefix_size`
pub struct KmerIterator<'a> {
    packed_slice: &'a [u32],
    k: usize,
    pos: usize,
    current_u64_window: u64,
    shift_inside_window: usize,
    filter_prefix_mask: u32,
}

impl<'a> KmerIterator<'a> {
    pub fn new(k: usize, packed_bytes: &'a [u128]) -> Self {
        let packed_slice = cast_slice(packed_bytes);
        Self {
            packed_slice,
            k,
            pos: 0,
            current_u64_window: get_u64_unaligned(packed_slice, 0),
            shift_inside_window: 0,
            filter_prefix_mask: u32::MAX >> (32 - 2 * k),
        }
    }
}

impl<'a> Iterator for KmerIterator<'a> {
    type Item = u32;
    fn next(&mut self) -> Option<Self::Item> {
        // If we reached the end of the current u64
        if self.shift_inside_window >= 32 {
            // End of the sequence
            if self.pos == self.packed_slice.len() - 2 {
                return None;
            }

            self.pos += 1;
            self.current_u64_window = get_u64_unaligned(self.packed_slice, self.pos);
            self.shift_inside_window = 0;
        }

        let current = ((self.current_u64_window >> self.shift_inside_window) as u32)
            & self.filter_prefix_mask;
        self.shift_inside_window += 2;
        Some(current)
    }
}

// Tests

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
            let kmer_it = KmerIterator::new(k, packed_bytes);
            ascii_bytes
                .windows(k)
                .zip(kmer_it)
                .for_each(|(ascii, kmer)| {
                    kmer_to_ascii(kmer, k, &mut buf);
                    assert_eq!(ascii, &buf[..k])
                });
        }
    }

    #[test]
    fn prefix_iterator() {
        let fasta: &[u8] = b">test_seq\nCCCTGAGTACGGAAAGCGCGAACGCAGATGCCCTATCGATACGTGGCAAGAGTGTTGTCCAAAGGGGCTACGCCCCTATTGAGTATTTACTATTGATTGTTAGATGTGAGTGCGTCTCAATCCTGCCGTTACTTGACCGTTTATGAGTTTATTATAGTCGTTAATATCTGGTCGAGACGGTGTAAAATACGCTATGCGACACCTGTCGTATATCAGAGAAAAGGGTCGATTCTCAATAAATATCGCCCTCTAAACCAGTTTAGGATGCTCTGGAGCCGAAGGATGGGTTCTTGCAGAATACATCACTTCTAGTAAGCGTCAGGCAAACGGCTTTAACCACCTTAGAAAGGGGCAATCACCCAAAGAATACAGTTGAGTAACGATTGTAAAAATAATGTAACAATGCATCAGTAGGAATCACCTTCACTTTCTTTGTATAGGAGTACGCACTCTTGTGGATACCCTCCGAACTACATACACGGTCCCAGTAACAGAGCT";
        for k in 1..=16 {
            check_kmer_iterator(fasta, k);
        }
    }
}
