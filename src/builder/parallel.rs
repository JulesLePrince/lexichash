use super::interleaved::InterleavedSlice32;
use super::single_thread::SingleThreadBuilder;
use crate::sketch::PartialSketch;
use crate::slice::SketchSlice32;
use crate::utils::*;
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
    pub fn build_with(&self, seq: &PackedDNA, res: &mut PartialSketch) {
        self.build_with_advanced::<true, true>(seq, res);
    }

    pub fn build_with_advanced<const PAR: bool, const PREFETCH: bool>(
        &self,
        seq: &PackedDNA,
        res: &mut PartialSketch,
    ) {
        let prefix_size = self.prefix_size;
        let suffix_size = self.k - self.prefix_size;
        let prefetch = PREFETCH && self.prefetch;
        let (packed_data, tail) = seq.bits();
        let packed_bytes = cast_slice::<u128, u8>(packed_data);
        let sketches = &mut res.sketches;
        res.num_kmers += seq.len() - (self.k - 1);

        if PAR {
            let missing_sketches = self.threads.saturating_sub(sketches.len());
            if missing_sketches > 0 {
                sketches.extend(
                    (0..missing_sketches).map(|_| InterleavedSlice32::from_masks(&self.masks)),
                );
            }

            let overlap = self.k.saturating_sub(1).div_ceil(4);
            let slices = overlapping_chunks(packed_bytes, self.threads, overlap);

            self.thread_pool.install(|| {
                slices
                    .into_par_iter()
                    .zip(sketches.par_iter_mut())
                    .for_each(|(&packed_bytes, res)| {
                        let builder = SingleThreadBuilder::new(prefix_size, suffix_size);
                        builder.build_with_dyn(packed_bytes, &mut res.0, prefetch);
                    });
            });
        } else {
            if sketches.is_empty() {
                sketches.push(InterleavedSlice32::from_masks(&self.masks));
            }
            let builder = SingleThreadBuilder::new(prefix_size, suffix_size);
            builder.build_with_dyn(packed_bytes, &mut sketches[0].0, prefetch);
        }

        // Tail processing
        let tail_nb_bases = seq.len() % 64;
        let packed_bytes_len = packed_bytes.len();
        let bytes_before_tail = unsafe {
            (packed_bytes[(packed_bytes_len - 8)..])
                .as_ptr()
                .cast::<u64>()
                .read_unaligned()
        };
        let mut first_window: u128 =
            (tail << (2 * (self.k - 1))) | ((bytes_before_tail >> (64 - 2 * (self.k - 1))) as u128);
        let mut second_window: u64 = (tail >> (128 - 2 * (self.k - 1))) as u64;

        let prefix_mask = u64::MAX >> (2 * (32 - self.prefix_size));
        let suffix_mask = u64::MAX >> (2 * (32 - self.k - self.prefix_size));

        for _ in 0..tail_nb_bases {
            // kmer processing
            let prefix = (first_window as u64) & prefix_mask;
            let suffix = (first_window as u64 >> (2 * self.prefix_size)) & suffix_mask;
            let s: u64 = (sketches[0].0[prefix as usize] >> 32) ^ suffix; // score of current kmer
            let best: u64 = u64::min(sketches[0].0[prefix as usize] >> 32, s);
            sketches[0].write_res(prefix as usize, best as u32);
            // rolling
            first_window = (first_window >> 2) | (((second_window & 0b11) as u128) << 126);
            second_window >>= 2;
        }
    }
}
