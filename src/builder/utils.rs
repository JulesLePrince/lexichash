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

pub fn packed_to_string(packed_bytes: u32, l: usize) -> String {
    let mut tmp = packed_bytes;
    let mut res = String::new();
    for _ in 0..l {
        let m: u32 = 0b11;
        let c = match tmp & m {
            0b00 => 'A',
            0b01 => 'C',
            0b10 => 'T',
            0b11 => 'G',
            _ => '_',
        };
        res.push(c);
        tmp >>= 2;
    }
    res
}
