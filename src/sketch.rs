use crate::builder::interleaved::InterleavedSlice32;
use crate::estimator::MutationRateEstimator;
use crate::slice::SketchSlice32;
use epserde::prelude::*;
use std::path::Path;
use wide::u32x8;

/// Accumulates the per-thread partial sketches and the total number of k-mers
/// seen across successive calls to [`SketchBuilder::build_with_advanced`](crate::SketchBuilder::build_with_advanced),
/// ready to be consolidated by [`PartialSketch::merge`].
pub struct PartialSketch {
    k: usize,
    prefix_size: usize,
    pub(crate) num_kmers: usize,
    pub(crate) sketches: Vec<InterleavedSlice32>,
}

impl PartialSketch {
    pub fn new(k: usize, prefix_size: usize) -> Self {
        Self {
            k,
            prefix_size,
            num_kmers: 0,
            sketches: Vec::new(),
        }
    }

    pub fn merge(&self) -> LexicSketch {
        let sketches = &self.sketches;
        let mut res = sketches[0].deinterleave();
        sketches.iter().skip(1).for_each(|sketch| {
            res.0
                .iter_mut()
                .zip(sketch.0.iter())
                .for_each(|(r, &c)| *r = (*r).min(c.1));
        });
        LexicSketch::new(self.k, self.prefix_size, res, self.num_kmers)
    }
}

#[derive(Epserde, Debug)]
pub struct LexicSketch {
    k: u8,
    prefix_size: u8,
    sketch_slice: SketchSlice32,
    num_kmers: usize,
}

impl LexicSketch {
    pub fn new(
        k: usize,
        prefix_size: usize,
        sketch_slice: SketchSlice32,
        num_kmers: usize,
    ) -> Self {
        Self {
            k: k as u8,
            prefix_size: prefix_size as u8,
            sketch_slice,
            num_kmers,
        }
    }

    pub fn get_k(&self) -> u8 {
        self.k
    }

    pub fn get_prefix_size(&self) -> u8 {
        self.prefix_size
    }

    pub fn get_sketch_slice(&self) -> &SketchSlice32 {
        &self.sketch_slice
    }

    pub fn serialize<P: AsRef<Path>>(&self, output_path: P) {
        unsafe { self.store(output_path).expect("Failed to lexichash data") }
    }

    pub fn deserialize<P: AsRef<Path>>(input_path: P) -> LexicSketch {
        // Fully allocates and copies the file into memory
        unsafe { <LexicSketch>::load_full(input_path).expect("Failed to load fully") }
    }

    pub fn average_match_size<'a>(&'a self, rhs: &'a Self) -> f64 {
        let prefix_size = self.prefix_size as usize;
        let suffix_size = self.k as usize - prefix_size;
        let padding = u32::BITS as usize / 2 - suffix_size;
        // padding prefix will be subtracted at the end
        let offset = prefix_size as f64 - padding as f64;
        self.sketch_slice
            .iter_leading_zeros(&rhs.sketch_slice)
            .fold(u32x8::ZERO, |u, v| u + (v >> 1))
            .reduce_add() as f64
            / (1 << (2 * prefix_size)) as f64
            + offset
    }

    #[inline(always)]
    pub fn get_divergence(&self, sk: &Self) -> f64 {
        self.get_divergence_with(sk, &self.make_estimator())
    }

    #[inline(always)]
    pub fn get_divergence_with(&self, sk: &Self, est: &MutationRateEstimator) -> f64 {
        Self::get_divergence_from_mean_with(self.average_match_size(sk), est)
    }

    #[inline(always)]
    pub fn get_divergence_from_mean(&self, mean: f64) -> f64 {
        Self::get_divergence_from_mean_with(mean, &self.make_estimator())
    }

    #[inline(always)]
    pub fn get_divergence_from_mean_with(mean: f64, est: &MutationRateEstimator) -> f64 {
        est.estimate_mut_rate::<2>(mean)
    }

    #[inline(always)]
    fn make_estimator(&self) -> MutationRateEstimator {
        MutationRateEstimator::new(self.k as usize, self.num_kmers)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use core::hint::black_box;

    #[test]
    #[ignore = "This is a benchmark, not a test"]
    fn bench_average_score() {
        let prefix_size = 6;
        let suffix_size = 16;
        let k = prefix_size + suffix_size;
        let rep = 1_000_000_000 >> (2 * prefix_size);

        let sketch1 = LexicSketch::new(
            k,
            prefix_size,
            SketchSlice32::random(prefix_size, suffix_size, 1),
            1_000_000,
        );
        let sketch2 = LexicSketch::new(
            k,
            prefix_size,
            SketchSlice32::random(prefix_size, suffix_size, 2),
            1_000_000,
        );

        let start = std::time::Instant::now();
        for _ in 0..rep {
            let res = sketch1.average_match_size(&sketch2);
            black_box(res);
        }
        eprintln!(
            "Computation of average score: {:.03} Gbp/s",
            start.elapsed().as_secs_f64().recip()
        );
    }
}
