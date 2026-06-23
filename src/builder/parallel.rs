use super::single_thread::SingleThreadBuilder;
use crate::LexicSketch;
use crate::slice::SketchSlice32;
use bytemuck::cast_slice;
use helicase::dna_format::PackedDNA;
use wide::u32x8;

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

        let (packed_data, tail) = seq.bits();
        let packed_bytes = cast_slice::<u128, u8>(packed_data);
        let single_thread_builder =
            SingleThreadBuilder::new(prefix_size, suffix_size, &self.masks.0);
        single_thread_builder.build_with::<true>(packed_bytes, &mut res[0].0);
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
