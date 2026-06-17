use super::iter::PrefixIter;
use super::utils::*;

/// Parameters for a sketch builder on one thread
pub struct SingleThreadBuilder<'a> {
    prefix_size: usize,
    suffix_size: usize,
    masks: &'a [u32],
}

impl<'a> SingleThreadBuilder<'a> {
    pub fn new(prefix_size: usize, suffix_size: usize, masks: &'a [u32]) -> Self {
        Self {
            prefix_size,
            suffix_size,
            masks,
        }
    }

    pub fn get_prefix_iterator(&self, packed_bytes: &'a [u128]) -> PrefixIter<'a> {
        PrefixIter::new(self.prefix_size, packed_bytes)
    }

    pub fn build(&self, packed_bytes: &'a [u128]) {
        for kmer in self.get_prefix_iterator(packed_bytes) {
            let seq_as_string = packed_to_string(kmer, self.prefix_size);
            println!("{seq_as_string}");
        }
    }
}
