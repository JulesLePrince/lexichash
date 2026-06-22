use super::iter::KmerIterator;

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

    #[inline(always)]
    pub fn get_prefix_iterator(&self, packed_bytes: &'a [u128]) -> KmerIterator<'a> {
        KmerIterator::new(self.prefix_size, packed_bytes)
    }

    #[inline(always)]
    pub fn get_suffix_iterator(&self, packed_bytes: &'a [u128]) -> impl Iterator<Item = u32> {
        KmerIterator::new(self.suffix_size, packed_bytes).skip(self.prefix_size)
    }

    /// Build while reusing an existing sketch, to avoid new allocation and merge on-the-fly
    pub fn build_with(&self, packed_bytes: &'a [u128], res: &'a mut [u32]) {
        let mut prefix_iterator = self.get_prefix_iterator(packed_bytes);
        let mut suffix_iterator = self.get_suffix_iterator(packed_bytes);
        while let Some(suffix) = suffix_iterator.next() {
            if let Some(prefix) = prefix_iterator.next() {
                let suffix_mask = self.masks[prefix as usize];
                // The min can be used without bias (to stick to lexichash definition it should be reverse min)
                res[prefix as usize] = u32::min(res[prefix as usize], suffix_mask ^ suffix);
            }
        }
    }
}
