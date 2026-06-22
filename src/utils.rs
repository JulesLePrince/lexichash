use std::io::{Write, stdout};

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

/// Prints a 2-bit packed k-mer to stdout without allocating.
pub fn print_kmer(kmer: u32, k: usize) {
    let mut buf = [0u8; 17];
    kmer_to_ascii(kmer, k, &mut buf);
    buf[k] = b'\n';
    stdout().write_all(&buf[..k + 1]).unwrap();
}
