use super::simd_iter::SimdKmerIterator;
use branches::{prefetch_read_data, prefetch_write_data};
use wide::u32x8;

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
    pub fn get_prefix_iterator(&self, packed_bytes: &'a [u8]) -> impl Iterator<Item = u32x8> {
        SimdKmerIterator::new(self.prefix_size, packed_bytes, 0)
    }

    #[inline(always)]
    pub fn get_suffix_iterator(&self, packed_bytes: &'a [u8]) -> impl Iterator<Item = u32x8> {
        SimdKmerIterator::new(self.suffix_size, packed_bytes, self.prefix_size)
    }

    #[inline(always)]
    pub fn build_with_dyn(&self, packed_bytes: &'a [u8], res: &'a mut [u32], prefetch: bool) {
        if prefetch {
            self.build_with::<true>(packed_bytes, res);
        } else {
            self.build_with::<false>(packed_bytes, res);
        }
    }

    pub fn build_with<const PREFETCH: bool>(&self, packed_bytes: &'a [u8], res: &'a mut [u32]) {
        let mut simd_prefix_iterator = self.get_prefix_iterator(packed_bytes);
        let mut simd_suffix_iterator = self.get_suffix_iterator(packed_bytes);
        // zip iterators together
        let simd_iter = simd_prefix_iterator
            .by_ref()
            .zip(simd_suffix_iterator.by_ref());

        // Previous vars initialized to have 0 effects on first iteration
        let mut prev_prefixes: [u32; 8] = [0; 8];
        let mut v_prev_suffix: u32x8 = u32x8::from([self.masks[0] ^ u32::MAX; 8]);

        for (v_prefix, v_suffix) in simd_iter {
            let prefixes = v_prefix.to_array();
            // ------ Prefetch ------
            if PREFETCH {
                prefixes.iter().map(|&p| p as usize).for_each(|p| {
                    prefetch_read_data::<_, 0>(&self.masks[p]);
                    prefetch_write_data::<_, 0>(&res[p]);
                });
            }

            // ------ Calculations ------
            let mut prev_masks = [0u32; 8];
            let mut prev_currents = [0u32; 8];

            for i in 0..8 {
                let prev_p = prev_prefixes[i] as usize;
                prev_masks[i] = self.masks[prev_p];
                prev_currents[i] = res[prev_p];
            }

            // Load into simd
            let v_prev_masks = u32x8::from(prev_masks);
            let v_prev_current = u32x8::from(prev_currents);

            // Operations
            let v_prev_result = v_prev_current.min(v_prev_masks ^ v_prev_suffix);
            let prev_results = v_prev_result.to_array();

            for i in 0..8 {
                res[prev_prefixes[i] as usize] = prev_results[i];
            }

            prev_prefixes = prefixes;
            v_prev_suffix = v_suffix;
        }
        // TODO : tail processing
    }
}
