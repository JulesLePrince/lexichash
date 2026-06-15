use crate::sketch::LexicSketch;
pub use helicase::dna_format::PackedDNA;

pub struct SketchBuilder {
    k: usize,
    prefix_size: usize,
    nb_threads: usize,
}

impl SketchBuilder {
    pub fn new(k: usize, prefix_size: usize, nb_threads: usize) -> Self {
        Self {
            k,
            prefix_size,
            nb_threads,
        }
    }

    pub fn build(&self, seq: &PackedDNA) -> LexicSketch {
        todo!();
    }
}
