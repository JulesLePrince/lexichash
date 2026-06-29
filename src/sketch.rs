use crate::slice::SketchSlice32;
use epserde::prelude::*;
use std::path::Path;
use wide::u32x8;

#[derive(Epserde, Debug)]
pub struct LexicSketch {
    k: u8,
    prefix_size: u8,
    sketch_slice: SketchSlice32,
}

impl LexicSketch {
    pub fn new(k: usize, prefix_size: usize, sketch_slice: SketchSlice32) -> Self {
        Self {
            k: k as u8,
            prefix_size: prefix_size as u8,
            sketch_slice,
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
        let offset = prefix_size as f64 - (u32::BITS / 2 - suffix_size as u32) as f64;
        self.sketch_slice
            .iter_leading_zeros(&rhs.sketch_slice)
            .fold(u32x8::ZERO, |u, v| u + (v >> 1))
            .reduce_add() as f64
            / (1 << (2 * prefix_size)) as f64
            + offset
    }

    pub fn get_score(&self, sk: &Self) -> f64 {
        todo!();
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
        );
        let sketch2 = LexicSketch::new(
            k,
            prefix_size,
            SketchSlice32::random(prefix_size, suffix_size, 2),
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
