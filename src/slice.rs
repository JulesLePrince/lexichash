use core::mem::transmute;
use epserde::Epserde;
use rand_xoshiro::Xoshiro256PlusPlus;
use rand_xoshiro::rand_core::{Rng, SeedableRng};
use wide::u32x8;

#[derive(Epserde, Debug, Clone)]
pub struct SketchSlice32(pub Vec<u32>);

impl SketchSlice32 {
    pub fn new(prefix_size: usize) -> Self {
        Self(vec![u32::MAX; 1 << (2 * prefix_size)])
    }

    #[allow(clippy::uninit_vec)]
    pub unsafe fn new_uninit(prefix_size: usize) -> Self {
        let len = 1 << (2 * prefix_size);
        let mut res: Vec<u32> = Vec::with_capacity(len);
        unsafe { res.set_len(len) };
        Self(res)
    }

    #[inline(always)]
    pub fn clear(&mut self) {
        self.0.fill(u32::MAX);
    }

    pub fn random(prefix_size: usize, suffix_size: usize, seed: u64) -> Self {
        let mut res = unsafe { Self::new_uninit(prefix_size) };
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        let bytes: &mut [u8] = bytemuck::cast_slice_mut(res.0.as_mut_slice());
        rng.fill_bytes(bytes);
        let shift = u32::BITS.saturating_sub(2 * suffix_size as u32);
        if shift > 0 {
            let (chunks, _) = res.0.as_chunks_mut::<8>();
            chunks.iter_mut().for_each(|chunk| {
                *chunk = (u32x8::new(*chunk) >> shift).to_array();
            });
        }
        res
    }

    /// Lazily yields 8 u32s at a time, without allocating.
    #[inline(always)]
    pub fn iter_chunks<'a>(&'a self) -> impl Iterator<Item = u32x8> + 'a {
        self.0.as_chunks::<8>().0.iter().copied().map(u32x8::new)
    }

    #[inline(always)]
    pub fn iter_leading_zeros<'a>(&'a self, rhs: &'a Self) -> impl Iterator<Item = u32x8> + 'a {
        debug_assert_eq!(self.0.len(), rhs.0.len());
        self.iter_chunks()
            .zip(rhs.iter_chunks())
            .map(|(u, v)| leading_zeros_u32x8(u ^ v))
    }

    /// Lazily yields the element-wise match size between `self` and `rhs`, 8 u32s at a time, without allocating.
    pub fn iter_matches<'a>(
        &'a self,
        rhs: &'a Self,
        prefix_size: usize,
        suffix_size: usize,
    ) -> impl Iterator<Item = u32x8> + 'a {
        let offset = (prefix_size as u32).wrapping_sub(u32::BITS / 2 - suffix_size as u32);
        let offsets = u32x8::new([offset; 8]);
        self.iter_leading_zeros(rhs)
            .map(move |lz| (lz >> 1) + offsets)
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

    use crate::utils::kmer_to_ascii;

    #[test]
    fn leading_zeros_matches_scalar() {
        let mut cases = vec![0u32, 1, 2, 3, u32::MAX, u32::MAX - 1, 0x80000000];
        cases.extend((0..32).map(|i| 1u32 << i));
        cases.extend((0..32).map(|i| u32::MAX >> i));
        cases.resize(cases.len().next_multiple_of(8), 0);

        for chunk in cases.chunks_exact(8) {
            let lanes: [u32; 8] = chunk.try_into().unwrap();
            let lz = leading_zeros_u32x8(u32x8::new(lanes)).to_array();
            let expected = lanes.map(u32::leading_zeros);
            assert_eq!(lz, expected, "input = {lanes:?}");
        }
    }

    fn check_iter_matches(prefix_size: usize, suffix_size: usize, seed: u64) {
        let sketch1 = SketchSlice32::random(prefix_size, suffix_size, seed + 1);
        let sketch2 = SketchSlice32::random(prefix_size, suffix_size, seed + 2);
        let (chunks1, _) = sketch1.0.as_chunks::<8>();
        let (chunks2, _) = sketch2.0.as_chunks::<8>();
        let mut buf1 = [0u8; 16];
        let mut buf2 = [0u8; 16];

        for (matches, (chunk1, chunk2)) in sketch1
            .iter_matches(&sketch2, prefix_size, suffix_size)
            .zip(chunks1.iter().zip(chunks2))
        {
            for ((&kmer1, &kmer2), &match_size) in chunk1.iter().zip(chunk2).zip(matches.as_array())
            {
                kmer_to_ascii(kmer1, suffix_size, &mut buf1);
                kmer_to_ascii(kmer2, suffix_size, &mut buf2);
                let common_suffix = buf1[..suffix_size]
                    .iter()
                    .rev()
                    .zip(buf2[..suffix_size].iter().rev())
                    .take_while(|(x, y)| x == y)
                    .count();
                let expected = (prefix_size + common_suffix) as u32;
                assert_eq!(match_size, expected);
            }
        }
    }

    #[test]
    fn test_iter_matches() {
        for prefix_size in 2..=7 {
            for suffix_size in 1..=16 {
                let seed = ((prefix_size * 16 + suffix_size - 1) * 2) as u64;
                check_iter_matches(prefix_size, suffix_size, seed);
            }
        }
    }
}
