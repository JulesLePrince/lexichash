use super::batched_iter::BatchedKmerIterator;
use branches::prefetch_write_data;
use wide::u32x8;

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
        const VECT: bool = cfg!(target_feature = "avx2");
        if prefetch {
            self.build_with::<true, VECT>(packed_bytes, mask_res);
        } else {
            self.build_with::<false, VECT>(packed_bytes, mask_res);
        }
    }

    pub fn build_with<const PREFETCH: bool, const VECT: bool>(
        &self,
        packed_bytes: &[u8],
        mask_res: &mut [u64],
    ) {
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
            Self::process_batch::<VECT>(mask_res, &prev_prefixes, &prev_suffixes);

            prev_prefixes = prefixes;
            prev_suffixes = suffixes;
        }
    }

    /// Gather the 8 buckets addressed by `prefixes`, update each running best with
    /// `best.min(mask ^ suffix)`, and scatter the results back.
    ///
    /// With `VECT` the per-lane `xor`/`min` is vectorized, otherwise it stays scalar.
    /// Vectorizing only paid off on AVX2.
    #[inline(always)]
    fn process_batch<const VECT: bool>(
        mask_res: &mut [u64],
        prefixes: &[u32; 8],
        suffixes: &[u32; 8],
    ) {
        let mut masks = [0u32; 8];
        let mut bests = [0u32; 8];
        for i in 0..8 {
            let c = unsafe { *mask_res.get_unchecked(prefixes[i] as usize) };
            masks[i] = c as u32;
            bests[i] = (c >> 32) as u32;
        }

        let results = if VECT {
            u32x8::new(bests)
                .min(u32x8::new(masks) ^ u32x8::new(*suffixes))
                .to_array()
        } else {
            let mut results = [0u32; 8];
            for i in 0..8 {
                results[i] = bests[i].min(masks[i] ^ suffixes[i]);
            }
            results
        };

        for i in 0..8 {
            unsafe {
                *mask_res.get_unchecked_mut(prefixes[i] as usize) =
                    ((results[i] as u64) << 32) | masks[i] as u64;
            }
        }
    }
}
