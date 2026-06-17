use super::utils::get_u64_unaligned;
use bytemuck::cast_slice;

/// Iterate `packed_bytes` over factors of size `prefix_size`
pub struct PrefixIter<'a> {
    packed_u32_bytes: &'a [u32],
    pos: usize,
    prefix_size: usize,
    current_u64_window: u64,
    remaining_inside_window: usize,
    shift_inside_window: usize,
    filter_prefix_mask: u32,
}

impl<'a> PrefixIter<'a> {
    pub fn new(prefix_size: usize, suffix_size: usize, packed_bytes: &'a [u128]) -> Self {
        let packed_u32_bytes = cast_slice(packed_bytes);
        Self {
            packed_u32_bytes,
            pos: 0,
            prefix_size,
            current_u64_window: get_u64_unaligned(packed_u32_bytes, suffix_size),
            remaining_inside_window: 17 - suffix_size,
            shift_inside_window: suffix_size,
            filter_prefix_mask: u32::MAX >> (32 - 2 * prefix_size),
        }
    }
}

impl<'a> Iterator for PrefixIter<'a> {
    type Item = u32;
    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining_inside_window == 0 {
            if self.pos == self.packed_u32_bytes.len() - 2 {
                return None;
            }

            self.pos += 1;
            self.current_u64_window = get_u64_unaligned(self.packed_u32_bytes, self.pos);
            self.remaining_inside_window = 17;
            self.shift_inside_window = 0;
        }

        let current = ((self.current_u64_window >> self.shift_inside_window) as u32)
            & self.filter_prefix_mask;
        self.shift_inside_window += 2;
        self.remaining_inside_window -= 1;
        Some(current)
    }
}

pub struct SuffixIter<'a> {
    packed_u32_bytes: &'a [u32],
    pos: usize,
    suffix_size: usize,
    current_u64_window: u64,
    remaining_inside_window: usize,
    shift_inside_window: usize,
    filter_suffix_mask: u32,
}


impl<'a> SuffixIter<'a> {
    pub fn new(prefix_size: usize, suffix_size: usize, packed_bytes: &'a [u128]) -> Self {
        let packed_u32_bytes = cast_slice(packed_bytes);
        Self {
            packed_u32_bytes,
            pos: 0,
            suffix_size,
            current_u64_window: get_u64_unaligned(packed_u32_bytes, 0),
            remaining_inside_window: 17,
            shift_inside_window: 0,
            filter_suffix_mask: u32::MAX >> (32 - 2 * suffix_size),
        }
    }
}

impl<'a> Iterator for SuffixIter<'a> {
    type Item = u32;
    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining_inside_window == 0 {
            if self.pos == self.packed_u32_bytes.len() - 2 {
                return None;
            }

            self.pos += 1;
            self.current_u64_window = get_u64_unaligned(self.packed_u32_bytes, self.pos);
            self.remaining_inside_window = 17;
            self.shift_inside_window = 0;
        }

        let current = ((self.current_u64_window >> self.shift_inside_window) as u32)
            & self.filter_suffix_mask;
        self.shift_inside_window += 2;
        self.remaining_inside_window -= 1;
        Some(current)
    }
}

// TODO add tests
