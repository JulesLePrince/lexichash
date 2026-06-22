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
            .fold(u32x8::ZERO, |u, v| u + v)
            .reduce_add() as f64
            / (1 << (2 * prefix_size + 1)) as f64
            + offset
    }

    pub fn get_score(&self, sk: &Self) -> f64 {
        todo!();
    }
}
