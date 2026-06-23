use super::batched_iter::BatchedKmerIterator;
use branches::prefetch_write_data;

/// Parameters for a sketch builder on one thread
pub struct SingleThreadBuilder {
    prefix_size: usize,
    suffix_size: usize,
}

impl SingleThreadBuilder {
    pub fn new(prefix_size: usize, suffix_size: usize) -> Self {
        Self {
            prefix_size,
            suffix_size,
        }
    }

    #[inline(always)]
    pub fn get_prefix_iterator<'a>(
        &self,
        packed_bytes: &'a [u8],
    ) -> impl Iterator<Item = [u32; 8]> + 'a {
        BatchedKmerIterator::new(self.prefix_size, packed_bytes, 0)
    }

    #[inline(always)]
    pub fn get_suffix_iterator<'a>(
        &self,
        packed_bytes: &'a [u8],
    ) -> impl Iterator<Item = [u32; 8]> + 'a {
        BatchedKmerIterator::new(self.suffix_size, packed_bytes, self.prefix_size)
    }

    #[inline(always)]
    pub fn build_with_dyn(&self, packed_bytes: &[u8], mask_res: &mut [u64], prefetch: bool) {
        if prefetch {
            self.build_with::<true>(packed_bytes, mask_res);
        } else {
            self.build_with::<false>(packed_bytes, mask_res);
        }
    }

    pub fn build_with<const PREFETCH: bool>(&self, packed_bytes: &[u8], mask_res: &mut [u64]) {
        assert_eq!(mask_res.len(), 1 << (2 * self.prefix_size));

        // Previous vars initialized to have 0 effects on first iteration
        let mut prev_prefixes: [u32; 8] = [0; 8];
        let mut prev_suffixes: [u32; 8] = [(mask_res[0] as u32) ^ u32::MAX; 8];

        for (prefixes, suffixes) in self
            .get_prefix_iterator(packed_bytes)
            .zip(self.get_suffix_iterator(packed_bytes))
        {
            // ------ Prefetch ------
            if PREFETCH {
                prefixes.iter().map(|&p| p as usize).for_each(|p| {
                    prefetch_write_data::<_, 0>(&mask_res[p]);
                });
            }

            // ------ Calculations ------
            let mut prev_masks = [0u32; 8];
            let mut prev_results = [0u32; 8];

            for i in 0..8 {
                let c = unsafe { *mask_res.get_unchecked(prev_prefixes[i] as usize) };
                let mask = c as u32;
                let best = (c >> 32) as u32;
                prev_masks[i] = mask;
                prev_results[i] = best.min(mask ^ prev_suffixes[i]);
            }

            for i in 0..8 {
                unsafe {
                    *mask_res.get_unchecked_mut(prev_prefixes[i] as usize) =
                        ((prev_results[i] as u64) << 32) | prev_masks[i] as u64;
                }
            }

            prev_prefixes = prefixes;
            prev_suffixes = suffixes;
        }
        // TODO : tail processing
    }
}
