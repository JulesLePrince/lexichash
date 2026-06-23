use smallvec::SmallVec;
use std::io::{Write, stdout};

/// Number of chunk references kept inline before spilling to the heap.
const INLINE_CHUNKS: usize = 16;

/// Splits `slice` into `num_chunks` contiguous subslices that each overlap the next one by `overlap` elements.
///
/// The slice is first evenly partitioned into `num_chunks` segments (the remainder is spread across the leading segments).
/// Every segment except the last is then extended by `overlap` elements into its successor, so consecutive subslices share `overlap` elements.
pub fn overlapping_chunks<T>(
    slice: &[T],
    num_chunks: usize,
    overlap: usize,
) -> SmallVec<[&[T]; INLINE_CHUNKS]> {
    let mut chunks = SmallVec::new();
    if num_chunks == 0 {
        return chunks;
    }

    let len = slice.len();
    let base = len / num_chunks;
    let rem = len % num_chunks;

    // Boundary of the non-overlapping partition that precedes chunk `i`,
    // distributing the `rem` leftover elements over the first `rem` chunks.
    let boundary = |i: usize| i * base + i.min(rem);

    for i in 0..num_chunks {
        let start = boundary(i);
        let end = (boundary(i + 1) + overlap).min(len);
        chunks.push(&slice[start..end]);
    }

    chunks
}

/// Get L1 data cache size in bytes.
///
/// Falls back to a conservative 32 KB when the platform cannot be queried.
pub fn l1_cache_bytes() -> usize {
    const FALLBACK: usize = 32 * 1024;

    #[cfg(target_os = "macos")]
    unsafe {
        let mut val: u64 = 0;
        let mut size = core::mem::size_of::<u64>();
        let rc = libc::sysctlbyname(
            c"hw.l1dcachesize".as_ptr(),
            &mut val as *mut u64 as *mut libc::c_void,
            &mut size,
            core::ptr::null_mut(),
            0,
        );
        if rc == 0 && val > 0 {
            return val as usize;
        }
    }

    #[cfg(target_os = "linux")]
    unsafe {
        let v = libc::sysconf(libc::_SC_LEVEL1_DCACHE_SIZE);
        if v > 0 {
            return v as usize;
        }
    }

    FALLBACK
}

pub fn get_u64_unaligned(sli: &[u32], i: usize) -> u64 {
    // TODO: avoid bound checks, this is slow

    // 1. Quick bounds check to ensure i + 1 exists (or use get() to avoid panics)
    if i + 1 >= sli.len() {
        panic!("Index out of bounds");
    }

    // 2. Get the raw pointer to the i-th element
    let ptr = unsafe { sli.as_ptr().add(i) as *const u64 };

    // 3. Read it as a u64
    unsafe { core::ptr::read_unaligned(ptr) }
}

/// Decodes the `k` low bases of a 2-bit packed k-mer into `buf` as ASCII.
pub fn kmer_to_ascii(kmer: u32, k: usize, buf: &mut [u8]) {
    const BASES: [u8; 4] = *b"ACTG";
    debug_assert!(k <= 16, "a u32 packs at most 16 bases");
    debug_assert!(k <= buf.len(), "`buf` must hold at least `k` bytes");

    let mut x = kmer;
    for slot in &mut buf[..k] {
        *slot = BASES[(x & 0b11) as usize];
        x >>= 2;
    }
}

pub fn kmer_to_string(kmer: u32, k: usize) -> String {
    let mut buf = [0u8; 32];
    kmer_to_ascii(kmer, k, &mut buf);
    String::from_utf8(buf[..k].to_vec()).expect("Invalid UTF-8 sequence")
}

/// Prints a 2-bit packed k-mer to stdout without allocating.
pub fn print_kmer(kmer: u32, k: usize) {
    let mut buf = [0u8; 17];
    kmer_to_ascii(kmer, k, &mut buf);
    buf[k] = b'\n';
    stdout().write_all(&buf[..k + 1]).unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chunks_cover_and_overlap() {
        let data: Vec<u32> = (0..10).collect();
        let chunks = overlapping_chunks(&data, 3, 2);

        // 10 elements over 3 chunks -> base 3, remainder 1 on the first chunk.
        // Boundaries: 0, 4, 7, 10. Each chunk extends 2 past its boundary.
        assert_eq!(chunks.len(), 3);
        assert_eq!(chunks[0], &[0, 1, 2, 3, 4, 5]); // [0..4) + 2
        assert_eq!(chunks[1], &[4, 5, 6, 7, 8]); // [4..7) + 2
        assert_eq!(chunks[2], &[7, 8, 9]); // [7..10), clamped

        // Consecutive chunks share exactly `overlap` elements where not clamped.
        assert_eq!(&chunks[0][chunks[0].len() - 2..], &chunks[1][..2]);
        assert_eq!(&chunks[1][chunks[1].len() - 2..], &chunks[2][..2]);
    }

    #[test]
    fn no_overlap_partitions_exactly() {
        let data: Vec<u32> = (0..9).collect();
        let chunks = overlapping_chunks(&data, 3, 0);
        assert_eq!(chunks.len(), 3);
        assert_eq!(chunks[0], &[0, 1, 2]);
        assert_eq!(chunks[1], &[3, 4, 5]);
        assert_eq!(chunks[2], &[6, 7, 8]);
    }

    #[test]
    fn zero_chunks_is_empty() {
        let data = [1, 2, 3];
        assert!(overlapping_chunks(&data, 0, 1).is_empty());
    }
}
