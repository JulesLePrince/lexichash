use super::simd_iter::SimdKmerIterator;
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
    pub fn get_prefix_iterator(&self, packed_bytes: &'a [u128]) -> impl Iterator<Item = u32x8> {
        SimdKmerIterator::new(self.prefix_size, packed_bytes, 0)
    }

    #[inline(always)]
    pub fn get_suffix_iterator(&self, packed_bytes: &'a [u128]) -> impl Iterator<Item = u32x8> {
        SimdKmerIterator::new(self.suffix_size, packed_bytes, self.prefix_size)
    }

    pub fn build_with(&self, packed_bytes: &'a [u128], res: &'a mut [u32]) {
        let mut simd_prefix_iterator = self.get_prefix_iterator(packed_bytes);
        let mut simd_suffix_iterator = self.get_suffix_iterator(packed_bytes);
        // zip iterators together
        let simd_iter = simd_prefix_iterator.by_ref().zip(simd_suffix_iterator.by_ref());
        for (v_prefix, v_suffix) in simd_iter {
            let prefixes = v_prefix.to_array();
            let mut masks = [0u32; 8];
            let mut currents = [0u32; 8];

            // Current values
            for i in 0..8 {
                let p = prefixes[i] as usize;
                masks[i] = self.masks[p];
                currents[i] = res[p];
            }

            // Load into simd
            let v_mask = u32x8::from(masks);
            let v_current = u32x8::from(currents);

            // Operations
            let v_result = v_current.min(v_mask ^ v_suffix);
            let results = v_result.to_array();

            for i in 0..8 {
                res[prefixes[i] as usize] = results[i];
            }
        }
        // TODO : tail processing
    }
}
