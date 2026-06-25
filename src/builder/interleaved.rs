use crate::slice::SketchSlice32;

/// A bucket table that interleaves the constant `mask` (low 32 bits)
/// and the running best value (high 32 bits) into a single `u64` per bucket.
///
/// Both are always read together in the build hot loop,
/// so keeping them in one `u64` lets each per-k-mer gather touch
/// a single cache line and issue one load instead of two.
/// Only the `best` half belongs to the sketch,
/// the mask half is discarded on [`deinterleave`](Self::deinterleave).
pub struct InterleavedSlice32(pub Vec<u64>);

impl InterleavedSlice32 {
    /// Build a fresh accumulator from a `masks` slice, with every `best` half initialized to `u32::MAX`.
    pub fn from_masks(masks: &SketchSlice32) -> Self {
        Self(
            masks
                .0
                .iter()
                .map(|&m| ((u32::MAX as u64) << 32) | m as u64)
                .collect(),
        )
    }

    #[inline(always)]
    pub fn write_res(&mut self, index: usize, res: u32) {
        let slot = unsafe { self.0.get_unchecked_mut(index) };
        *slot = ((res as u64) << 32) | (*slot & (u32::MAX as u64));
    }

    /// Drop the masks and return the `best`/`res` half as a `SketchSlice32`.
    pub fn deinterleave(&self) -> SketchSlice32 {
        SketchSlice32(self.0.iter().map(|&c| (c >> 32) as u32).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_masks_initializes_best_to_max() {
        let masks = SketchSlice32(vec![1, 2, 3, 4]);
        let interleaved = InterleavedSlice32::from_masks(&masks);
        // low halves hold the masks, high halves are all u32::MAX
        assert_eq!(
            interleaved.0,
            vec![
                (u32::MAX as u64) << 32 | 1,
                (u32::MAX as u64) << 32 | 2,
                (u32::MAX as u64) << 32 | 3,
                (u32::MAX as u64) << 32 | 4,
            ]
        );
        assert_eq!(interleaved.deinterleave().0, vec![u32::MAX; 4]);
    }

    #[test]
    fn deinterleave_extracts_high_halves() {
        let interleaved = InterleavedSlice32(vec![
            (10u64 << 32) | 1,
            (20u64 << 32) | 2,
            (30u64 << 32) | 3,
        ]);
        assert_eq!(interleaved.deinterleave().0, vec![10, 20, 30]);
    }
}
