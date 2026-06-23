use super::interleaved::InterleavedSlice32;
use super::single_thread::SingleThreadBuilder;
use crate::LexicSketch;
use crate::slice::SketchSlice32;
use crate::utils::{l1_cache_bytes, overlapping_chunks};
use bytemuck::cast_slice;
use helicase::dna_format::PackedDNA;
use rayon::prelude::*;
use rayon::{ThreadPool, ThreadPoolBuilder};

/// The Sketch Builder Parameters
pub struct SketchBuilder {
    k: usize,
    prefix_size: usize,
    threads: usize,
    masks: SketchSlice32,
    thread_pool: ThreadPool,
    prefetch: bool,
}

impl SketchBuilder {
    pub fn new(k: usize, prefix_size: usize, threads: usize) -> Self {
        let suffix_size = k - prefix_size;
        let masks = SketchSlice32::random(prefix_size, suffix_size, 101010);
        Self::new_with_masks(k, prefix_size, threads, masks)
    }

    pub fn new_with_masks(
        k: usize,
        prefix_size: usize,
        threads: usize,
        masks: SketchSlice32,
    ) -> Self {
        let thread_pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("Failed to build thread pool");
        let table_bytes = (1usize << (2 * prefix_size)) * core::mem::size_of::<u32>();
        let prefetch = 2 * table_bytes > l1_cache_bytes();
        Self {
            k,
            prefix_size,
            threads,
            masks,
            thread_pool,
            prefetch,
        }
    }

    #[inline(always)]
    pub fn build_with(&self, seq: &PackedDNA, res: &mut Vec<InterleavedSlice32>) {
        self.build_with_advanced::<true, true>(seq, res);
    }

    pub fn build_with_advanced<const PAR: bool, const PREFETCH: bool>(
        &self,
        seq: &PackedDNA,
        res: &mut Vec<InterleavedSlice32>,
    ) {
        let prefix_size = self.prefix_size;
        let suffix_size = self.k - self.prefix_size;
        let prefetch = PREFETCH && self.prefetch;
        let (packed_data, tail) = seq.bits();
        let packed_bytes = cast_slice::<u128, u8>(packed_data);

        if PAR {
            let missing_sketches = self.threads.saturating_sub(res.len());
            if missing_sketches > 0 {
                res.extend(
                    (0..missing_sketches).map(|_| InterleavedSlice32::from_masks(&self.masks)),
                );
            }

            let overlap = self.k.saturating_sub(1).div_ceil(4);
            let slices = overlapping_chunks(packed_bytes, self.threads, overlap);

            self.thread_pool.install(|| {
                slices
                    .into_par_iter()
                    .zip(res.par_iter_mut())
                    .for_each(|(&packed_bytes, res)| {
                        let builder = SingleThreadBuilder::new(prefix_size, suffix_size);
                        builder.build_with_dyn(packed_bytes, &mut res.0, prefetch);
                    });
            });
        } else {
            if res.is_empty() {
                res.push(InterleavedSlice32::from_masks(&self.masks));
            }
            let builder = SingleThreadBuilder::new(prefix_size, suffix_size);
            builder.build_with_dyn(packed_bytes, &mut res[0].0, prefetch);
        }

        // TODO tail processing
    }

    pub fn merge_sketches(&self, sketches: &[InterleavedSlice32]) -> LexicSketch {
        let mut res = sketches[0].deinterleave();
        sketches.iter().skip(1).for_each(|sketch| {
            res.0
                .iter_mut()
                .zip(sketch.0.iter())
                .for_each(|(r, &c)| *r = (*r).min((c >> 32) as u32));
        });
        LexicSketch::new(self.k, self.prefix_size, res)
    }
}
