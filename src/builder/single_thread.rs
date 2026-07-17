use super::batched_iter::BatchedKmerIterator;
use crate::builder::interleaved::InterleavedSlice32;
use branches::prefetch_write_data;
use wide::u32x8;

/// Parameters for a sketch builder on one thread
pub struct SingleThreadBuilder {
    prefix_size: usize,
    suffix_size: usize,
}

impl SingleThreadBuilder {
    pub const fn new(prefix_size: usize, suffix_size: usize) -> Self {
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
    pub fn process_seq(
        &self,
        packed_bytes: &[u8],
        mask_res: &mut InterleavedSlice32,
        prefetch: bool,
    ) {
        const VECT: bool = cfg!(any(target_feature = "avx2", target_feature = "neon"));
        if prefetch {
            self.process_seq_advanced::<true, VECT>(packed_bytes, mask_res);
        } else {
            self.process_seq_advanced::<false, VECT>(packed_bytes, mask_res);
        }
    }

    pub fn process_seq_advanced<const PREFETCH: bool, const VECT: bool>(
        &self,
        packed_bytes: &[u8],
        mask_res: &mut InterleavedSlice32,
    ) {
        assert_eq!(mask_res.len(), 1 << (2 * self.prefix_size));

        // Previous vars initialized to have 0 effects on first iteration
        let mut prev_prefixes: [u32; 8] = [0; 8];
        let mut prev_suffixes: [u32; 8] = [mask_res.get_mask(0) ^ u32::MAX; 8];

        for (prefixes, suffixes) in self
            .get_prefix_iterator(packed_bytes)
            .zip(self.get_suffix_iterator(packed_bytes))
        {
            // ------ Prefetch ------
            if PREFETCH {
                prefixes.iter().map(|&p| p as usize).for_each(|p| {
                    prefetch_write_data::<_, 0>(unsafe { mask_res.0.as_ptr().add(p) });
                });
            }

            // ------ Calculations ------
            Self::process_chunk::<VECT>(mask_res, &prev_prefixes, &prev_suffixes);

            prev_prefixes = prefixes;
            prev_suffixes = suffixes;
        }
    }

    /// Gather the 8 buckets addressed by `prefixes`, update each running best with
    /// `best.min(mask ^ suffix)`, and scatter the results back.
    ///
    /// With `VECT` the per-lane `xor`/`min` is vectorized, otherwise it stays scalar.
    #[inline(always)]
    fn process_chunk<const VECT: bool>(
        mask_res: &mut InterleavedSlice32,
        prefixes: &[u32; 8],
        suffixes: &[u32; 8],
    ) {
        let mut masks = [0u32; 8];
        let mut bests = [0u32; 8];
        for ((mask, best), &prefix) in masks.iter_mut().zip(bests.iter_mut()).zip(prefixes) {
            (*mask, *best) = mask_res.get_mask_res(prefix as usize);
        }

        let results = if VECT {
            u32x8::new(bests)
                .min(u32x8::new(masks) ^ u32x8::new(*suffixes))
                .to_array()
        } else {
            let mut results = [0u32; 8];
            for (result, (&best, (&mask, &suffix))) in results
                .iter_mut()
                .zip(bests.iter().zip(masks.iter().zip(suffixes)))
            {
                *result = best.min(mask ^ suffix);
            }
            results
        };

        for (&result, &prefix) in results.iter().zip(prefixes) {
            mask_res.set_res(prefix as usize, result);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::slice::SketchSlice32;
    use core::hint::black_box;
    use rand_xoshiro::Xoshiro256PlusPlus;
    use rand_xoshiro::rand_core::{Rng, SeedableRng};

    fn random_packed_bytes(num_bases: usize, seed: u64) -> Vec<u8> {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        let mut bytes = vec![0u8; num_bases / 4];
        rng.fill_bytes(&mut bytes);
        bytes
    }

    fn bench_variant<const VECT: bool>(
        name: &str,
        packed_bytes: &[u8],
        builder: &SingleThreadBuilder,
        masks: &SketchSlice32,
        rep: usize,
    ) {
        let bases = packed_bytes.len() * 4;
        let start = std::time::Instant::now();
        for _ in 0..rep {
            let mut mask_res = InterleavedSlice32::from_masks(masks);
            builder.process_seq_advanced::<false, VECT>(black_box(packed_bytes), &mut mask_res);
            black_box(&mask_res);
        }
        let elapsed = start.elapsed().as_secs_f64();
        let gbp_s = (bases * rep) as f64 / elapsed / 1e9;
        eprintln!("{name}: {gbp_s:.03} Gbp/s ({elapsed:.03}s total, {rep} reps)");
    }

    #[test]
    #[ignore = "This is a benchmark, not a test"]
    fn bench_vect_neon() {
        let prefix_size = 6;
        let suffix_size = 16;
        let num_bases = 10_000_000;
        let rep = 30;

        let packed_bytes = random_packed_bytes(num_bases, 1);
        let masks = SketchSlice32::random(prefix_size, suffix_size, 2);
        let builder = SingleThreadBuilder::new(prefix_size, suffix_size);

        bench_variant::<true>("SIMD (VECT=true)", &packed_bytes, &builder, &masks, rep);
        bench_variant::<false>("scalar (VECT=false)", &packed_bytes, &builder, &masks, rep);
    }
}
