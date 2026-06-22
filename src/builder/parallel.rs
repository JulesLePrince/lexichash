use super::single_thread::SingleThreadBuilder;
use crate::LexicSketch;
use crate::slice::SketchSlice32;
use helicase::dna_format::PackedDNA;

/// The Sketch Builder Parameters
pub struct SketchBuilder {
    k: usize,
    prefix_size: usize,
    threads: usize,
    masks: SketchSlice32,
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
        Self {
            k,
            prefix_size,
            threads,
            masks,
        }
    }

    pub fn build_with(&self, seq: &PackedDNA, res: &mut Vec<SketchSlice32>) {
        let prefix_size = self.prefix_size;
        let suffix_size = self.k - self.prefix_size;
        let missing_sketches = self.threads.saturating_sub(res.len());
        if missing_sketches > 0 {
            res.extend((0..missing_sketches).map(|_| SketchSlice32::new(prefix_size)));
        }

        let (packed_bytes, _) = seq.bits();
        let single_thread_builder =
            SingleThreadBuilder::new(prefix_size, suffix_size, &self.masks.0);
        single_thread_builder.build_with(packed_bytes, &mut res[0].0);
    }

    pub fn merge_sketches(&self, sketches: &[SketchSlice32]) -> LexicSketch {
        todo!()
    }
}
