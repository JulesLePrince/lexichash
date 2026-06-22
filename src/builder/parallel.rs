use super::single_thread::SingleThreadBuilder;
use crate::LexicSketch;
use crate::slice::SketchSlice32;
use helicase::dna_format::PackedDNA;

/// The Sketch Builder Parameters
pub struct SketchBuilder {
    k: usize,
    prefix_size: usize,
    threads: usize,
}

impl SketchBuilder {
    pub fn new(k: usize, prefix_size: usize, threads: usize) -> Self {
        Self {
            k,
            prefix_size,
            threads,
        }
    }

    pub fn build(&self, seq: &PackedDNA) -> LexicSketch {
        let suffix_size = self.k - self.prefix_size;
        let masks = SketchSlice32::random(self.prefix_size, 2, 101010);
        let (packed_bytes, _) = seq.bits();
        let single_thread_builder =
            SingleThreadBuilder::new(self.prefix_size, suffix_size, &masks.0);
        let fingerprint = single_thread_builder.build(packed_bytes);
        LexicSketch::new(self.k, self.prefix_size, SketchSlice32(fingerprint))
    }
}
