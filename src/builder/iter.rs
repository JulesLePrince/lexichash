use super::utils::get_u64_unaligned;
use bytemuck::cast_slice;

/// Iterate `packed_bytes` over factors of size `prefix_size`
pub struct KmerIterator<'a> {
    packed_u32_bytes: &'a [u32],
    pos: usize,
    size: usize,
    current_u64_window: u64,
    shift_inside_window: usize,
    filter_prefix_mask: u32,
}

impl<'a> KmerIterator<'a> {
    pub fn new(size: usize, packed_bytes: &'a [u128]) -> Self {
        let packed_u32_bytes = cast_slice(packed_bytes);
        Self {
            packed_u32_bytes: packed_u32_bytes,
            pos: 0,
            size: size,
            current_u64_window: get_u64_unaligned(packed_u32_bytes, 0),
            shift_inside_window: 0,
            filter_prefix_mask: u32::MAX >> (32 - 2 * size),
        }
    }
}

impl<'a> Iterator for KmerIterator<'a> {
    type Item = u32;
    fn next(&mut self) -> Option<Self::Item> {
        // If we reached the end of the current u64
        if self.shift_inside_window >= 32 {
            // End of the sequence
            if self.pos == self.packed_u32_bytes.len() - 2 {
                return None;
            }

            self.pos += 1;
            self.current_u64_window = get_u64_unaligned(self.packed_u32_bytes, self.pos);
            self.shift_inside_window = 0;
        }

        let current = ((self.current_u64_window >> self.shift_inside_window) as u32)
            & self.filter_prefix_mask;
        self.shift_inside_window += 2;
        Some(current)
    }
}


// TODO add tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prefix_iterator() {
        // Generate random_seq
        let seq: &[u32] = &[0xDEADBEEF; 1024];

    }
}
