use crate::slice::SketchSlice32;
use std::{path::Path, slice::SliceIndex};

use epserde::prelude::*;

#[derive(Epserde, Debug)]
pub struct LexicSketch {
    k: u8,
    prefix_size: u8,
    fingerprint: SketchSlice32,
}

impl LexicSketch {
    pub fn new(k: u8, prefix_size: u8, fingerprint: SketchSlice32) -> Self {
        return Self { k, prefix_size, fingerprint }
    }

    pub fn get_k(&self) -> u8 {
        self.k
    }

    pub fn get_prefix_size(&self) -> u8 {
        self.prefix_size
    }

    pub fn get_fingerprint(&self) -> &SketchSlice32 {
        &self.fingerprint
    }

    pub fn serialize(&self, output_path: impl AsRef<Path>) -> () {
        unsafe {
            self.store(output_path).expect("Failed to lexichash data")
        }
    }

    pub fn deserialize(input_path: impl AsRef<Path>) -> LexicSketch {
        // Fully allocates and copies the file into memory
        unsafe {
            let full_data = <LexicSketch>::load_full(input_path).expect("Failed to load fully");
            return full_data;
        }
    }

    pub fn get_score(&self, sk: &Self) -> f64 {
        todo!();
    }
}
