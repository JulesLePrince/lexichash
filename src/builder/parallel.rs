use super::single_thread::SingleThreadBuilder;
use crate::LexicSketch;
use crate::slice::SketchSlice32;
use crate::utils::overlapping_chunks;
use bytemuck::cast_slice;
use helicase::dna_format::PackedDNA;
use rayon::prelude::*;
use rayon::{ThreadPool, ThreadPoolBuilder};
use wide::u32x8;

/// The Sketch Builder Parameters
pub struct SketchBuilder {
    k: usize,
    prefix_size: usize,
    threads: usize,
    masks: SketchSlice32,
    thread_pool: ThreadPool,
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
        Self {
            k,
            prefix_size,
            threads,
            masks,
            thread_pool,
        }
    }

    pub fn build_with_advanced<const PAR: bool, const PREFETCH: bool>(
        &self,
        seq: &PackedDNA,
        res: &mut Vec<SketchSlice32>,
    ) {
        let prefix_size = self.prefix_size;
        let suffix_size = self.k - self.prefix_size;
        let (packed_data, tail) = seq.bits();
        let packed_bytes = cast_slice::<u128, u8>(packed_data);

        if PAR {
            let missing_sketches = self.threads.saturating_sub(res.len());
            if missing_sketches > 0 {
                res.extend((0..missing_sketches).map(|_| SketchSlice32::new(prefix_size)));
            }

            let overlap = self.k.saturating_sub(1).div_ceil(4);
            let slices = overlapping_chunks(packed_bytes, self.threads, overlap);

            self.thread_pool.install(|| {
                slices
                    .into_par_iter()
                    .zip(res.par_iter_mut())
                    .for_each(|(&packed_bytes, res)| {
                        let builder =
                            SingleThreadBuilder::new(prefix_size, suffix_size, &self.masks.0);
                        builder.build_with::<PREFETCH>(packed_bytes, &mut res.0);
                    });
            });
        } else {
            if res.is_empty() {
                res.push(SketchSlice32::new(prefix_size));
            }
            let builder = SingleThreadBuilder::new(prefix_size, suffix_size, &self.masks.0);
            builder.build_with::<PREFETCH>(packed_bytes, &mut res[0].0);
        }

        // TODO tail processing
    }

    pub fn build_with(&self, seq: &PackedDNA, res: &mut Vec<SketchSlice32>) {
        let prefix_size = self.prefix_size;
        let suffix_size = self.k - self.prefix_size;
        let missing_sketches = self.threads.saturating_sub(res.len());
        if missing_sketches > 0 {
            res.extend((0..missing_sketches).map(|_| SketchSlice32::new(prefix_size)));
        }

        let (packed_data, tail) = seq.bits();
        let packed_bytes = cast_slice::<u128, u8>(packed_data);
        let overlap = self.k.saturating_sub(1).div_ceil(4);
        let slices = overlapping_chunks(packed_bytes, self.threads, overlap);

        self.thread_pool.install(|| {
            slices
                .into_par_iter()
                .zip(res.par_iter_mut())
                .for_each(|(&packed_bytes, res)| {
                    let builder = SingleThreadBuilder::new(prefix_size, suffix_size, &self.masks.0);
                    builder.build_with::<true>(packed_bytes, &mut res.0);
                });
        });

        // TODO tail processing
    }

    pub fn merge_sketches(&self, sketches: &[SketchSlice32]) -> LexicSketch {
        let mut res = sketches[0].clone();
        let (chunks, _) = res.0.as_chunks_mut::<8>();
        sketches.iter().skip(1).for_each(|sketch| {
            chunks
                .iter_mut()
                .zip(sketch.iter_chunks())
                .for_each(|(chunk, v)| {
                    *chunk = (u32x8::new(*chunk).min(v)).to_array();
                })
        });
        LexicSketch::new(self.k, self.prefix_size, res)
    }
}
