use crate::slice::SketchSlice32;

/// A bucket table that interleaves the constant `mask` (low 32 bits)
/// and the running best value (high 32 bits) into a single `u64` per bucket.
///
/// Both are always read together in the build hot loop,
/// so keeping them in one `u64` lets each per-k-mer gather touch
/// a single cache line and issue one load instead of two.
/// Only the `best` half belongs to the sketch,
/// the mask half is discarded on [`deinterleave`](Self::deinterleave).
pub struct InterleavedSlice32(pub Vec<(u32, u32)>);

impl InterleavedSlice32 {
    /// Build a fresh accumulator from a `masks` slice, with every `best` half initialized to `u32::MAX`.
    pub fn from_masks(masks: &SketchSlice32) -> Self {
        Self(masks.0.iter().map(|&m| (m, u32::MAX)).collect())
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    #[allow(unused)]
    #[inline(always)]
    pub fn get_mask_res(&self, index: usize) -> (u32, u32) {
        unsafe { *self.0.get_unchecked(index) }
    }

    #[inline(always)]
    pub fn get_mask(&self, index: usize) -> u32 {
        unsafe { self.0.get_unchecked(index).0 }
    }

    #[inline(always)]
    pub fn get_res(&self, index: usize) -> u32 {
        unsafe { self.0.get_unchecked(index).1 }
    }

    #[inline(always)]
    pub fn set_res(&mut self, index: usize, res: u32) {
        unsafe { self.0.get_unchecked_mut(index).1 = res }
    }

    /// Drop the masks and return the `best`/`res` half as a `SketchSlice32`.
    pub fn deinterleave(&self) -> SketchSlice32 {
        SketchSlice32(self.0.iter().map(|&c| c.1).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_masks_initializes_best_to_max() {
        let masks = SketchSlice32(vec![1, 2, 3, 4]);
        let interleaved = InterleavedSlice32::from_masks(&masks);
        assert_eq!(
            interleaved.0,
            vec![(1, u32::MAX), (2, u32::MAX), (3, u32::MAX), (4, u32::MAX)]
        );
        assert_eq!(interleaved.deinterleave().0, vec![u32::MAX; 4]);
    }

    #[test]
    fn deinterleave_extracts_best_halves() {
        let interleaved = InterleavedSlice32(vec![(1, 10), (2, 20), (3, 30)]);
        assert_eq!(interleaved.deinterleave().0, vec![10, 20, 30]);
    }
}
