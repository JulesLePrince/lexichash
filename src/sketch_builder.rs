use crate::sketch::LexicSketch;
use helicase::dna_format::PackedDNA;

struct SketchBuilder {
    k: usize,
    prefix_size: usize,
    nb_threads: usize,
}

impl SketchBuilder  {
    fn init(k: usize, prefix_size: usize, nb_threads: usize) -> Self {
        Self {
            k,
            prefix_size,
            nb_threads,
        }
    }

    fn build(&self, seq: &PackedDNA) -> LexicSketch {
        todo!();
    }
}
