use crate::builder::utils;

use super::iter::KmerIterator;
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

    pub fn get_prefix_iterator(&self, packed_bytes: &'a [u128]) -> KmerIterator<'a> {
        KmerIterator::new(self.prefix_size, packed_bytes)
    }

    pub fn get_suffix_iterator(&self, packed_bytes: &'a [u128]) -> KmerIterator<'a> {
        KmerIterator::new(self.suffix_size, packed_bytes)
    }


    pub fn build(&self, packed_bytes: &'a [u128]) {
        let mut prefix_iterator = self.get_prefix_iterator(packed_bytes);
        let mut suffix_iterator = self.get_suffix_iterator(packed_bytes);
        suffix_iterator.nth(self.prefix_size-1);
        while let Some(suffix) = suffix_iterator.next() {
            if let Some(prefix) = prefix_iterator.next() {
                println!("{}{}", utils::packed_to_string(prefix, self.prefix_size), utils::packed_to_string(suffix, self.suffix_size))
            }
        }
    }
}
