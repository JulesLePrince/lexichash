use core::mem::transmute;
use rand_xoshiro::Xoshiro256PlusPlus;
use rand_xoshiro::rand_core::{Rng, SeedableRng};
use wide::u32x8;
use epserde::Epserde;

#[derive(Epserde, Debug)]
pub struct SketchSlice32(pub Vec<u32>);

impl SketchSlice32 {
    pub fn new(prefix_size: usize) -> Self {
        Self(vec![u32::MAX; 1 << (2 * prefix_size)])
    }

    #[allow(clippy::uninit_vec)]
    pub fn random(prefix_size: usize, suffix_size: usize, seed: u64) -> Self {
        let len = 1 << (2 * prefix_size);
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        let mut res: Vec<u32> = Vec::with_capacity(len);
        unsafe { res.set_len(len) };
        let bytes: &mut [u8] = bytemuck::cast_slice_mut(res.as_mut_slice());
        rng.fill_bytes(bytes);
        let num_bits = (2 * suffix_size) as u32;
        let num_bits_mask = u32::MAX >> (32 - 2*(num_bits));
        for val in res.iter_mut() {
            *val = *val & num_bits_mask;
        }
        Self(res)
    }

    /// Lazily yields the element-wise XOR of `self` and `rhs`, 8 u32s at a time, without allocating.
    pub fn xor_chunks<'a>(&'a self, rhs: &'a SketchSlice32) -> impl Iterator<Item = u32x8> + 'a {
        assert_eq!(self.0.len(), rhs.0.len());
        let (lhs, _) = self.0.as_chunks::<8>();
        let (rhs, _) = rhs.0.as_chunks::<8>();
        lhs.iter()
            .zip(rhs)
            .map(|(u, v)| u32x8::new(*u) ^ u32x8::new(*v))
    }
}

/// Per-lane count of leading zero bits.
#[inline(always)]
#[cfg(target_feature = "neon")]
pub fn leading_zeros_u32x8(x: u32x8) -> u32x8 {
    use core::arch::aarch64::{uint32x4_t, vclzq_u32};
    unsafe {
        let [lo, hi]: [uint32x4_t; 2] = transmute(x);
        transmute([vclzq_u32(lo), vclzq_u32(hi)])
    }
}

/// Per-lane count of leading zero bits.
#[inline(always)]
#[cfg(target_feature = "avx2")]
pub fn leading_zeros_u32x8(x: u32x8) -> u32x8 {
    use core::arch::x86_64::*;
    unsafe {
        let mut v: __m256i = transmute(x);
        v = _mm256_or_si256(v, _mm256_srli_epi32::<1>(v));
        v = _mm256_or_si256(v, _mm256_srli_epi32::<2>(v));
        v = _mm256_or_si256(v, _mm256_srli_epi32::<4>(v));
        v = _mm256_or_si256(v, _mm256_srli_epi32::<8>(v));
        v = _mm256_or_si256(v, _mm256_srli_epi32::<16>(v));
        let lut = _mm256_setr_epi8(
            0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, //
            0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
        );
        let low_mask = _mm256_set1_epi8(0x0f);
        let lo = _mm256_and_si256(v, low_mask);
        let hi = _mm256_and_si256(_mm256_srli_epi16::<4>(v), low_mask);
        // byte-level popcount
        let bytes = _mm256_add_epi8(_mm256_shuffle_epi8(lut, lo), _mm256_shuffle_epi8(lut, hi));
        // sum the 4 bytes within each 32-bit lane via `(b * 0x01010101) >> 24`
        let sum = _mm256_srli_epi32::<24>(_mm256_mullo_epi32(bytes, _mm256_set1_epi32(0x01010101)));
        transmute(_mm256_sub_epi32(_mm256_set1_epi32(32), sum))
    }
}

/// Per-lane count of leading zero bits.
#[inline(always)]
#[cfg(not(any(target_feature = "neon", target_feature = "avx2")))]
pub fn leading_zeros_u32x8(x: u32x8) -> u32x8 {
    u32x8::new(x.to_array().map(u32::leading_zeros))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn leading_zeros_matches_scalar() {
        let mut cases = vec![0u32, 1, 2, 3, u32::MAX, u32::MAX - 1, 0x80000000];
        cases.extend((0..32).map(|i| 1u32 << i));
        cases.extend((0..32).map(|i| u32::MAX >> i));
        cases.resize(cases.len().next_multiple_of(8), 0);

        for chunk in cases.chunks_exact(8) {
            let lanes: [u32; 8] = chunk.try_into().unwrap();
            let got = leading_zeros_u32x8(u32x8::new(lanes)).to_array();
            let want = lanes.map(u32::leading_zeros);
            assert_eq!(got, want, "input = {lanes:?}");
        }
    }
}
