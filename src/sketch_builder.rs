use std::io::Cursor;

use crate::sketch::LexicSketch;
use helicase::dna_format::PackedDNA;
use fastrand;


// ---------------------
//  Struct definitions
// ---------------------

/// The Sketch Builder Parameters
pub struct SketchBuilder {
    k: usize,
    prefix_size: usize,
    nb_threads: usize,
}

/// Parameters for a sketch builder on one thread
struct SingleThreadedBuilder{
    masks: Vec<u32>,
    k: usize,
    prefix_size: usize,
}

/// Iterate `packed_bytes` over factors of size `prefix_size`
struct PrefixIter<'a> {
    packed_u32_bytes: &'a [u32],
    pos: usize,
    prefix_size: usize,
    current_u64_window: u64,
    remaining_inside_window: usize,
    shift_inside_window: usize,
    filter_prefix_mask: u32,
}

// ---------------------
//        Utils
// ---------------------

fn create_masks(skb: &SketchBuilder) -> Vec<u32> {
    let nb_masks = usize::pow(4, skb.prefix_size as u32);
    let suffix_size = skb.k - skb.prefix_size;
    let mut res: Vec<u32> = vec![0; nb_masks];
    for i in 0..nb_masks {
        // We shift to have zeros, to prevent xor bugs
        res[i] = fastrand::u32(..) >> (32-2*(suffix_size));
    }
    return res
}

fn transform_to_u32_slice_unsafe(source: &[u128]) -> &[u32] {
    let ptr = source.as_ptr() as *const u32;
    let len = source.len() * 4; // 1 u128 becomes 4 u32s
    unsafe {
        std::slice::from_raw_parts(ptr, len)
    }
}

pub fn get_u64_unaligned(sli: &[u32], i: usize) -> u64 {
    // 1. Quick bounds check to ensure i + 1 exists (or use get() to avoid panics)
    if i + 1 >= sli.len() {
        panic!("Index out of bounds");
    }

    // 2. Get the raw pointer to the i-th element
    let ptr = unsafe { sli.as_ptr().add(i) as *const u64 };

    // 3. Read it as a u64
    unsafe { core::ptr::read_unaligned(ptr) }
}

fn packed_to_string(packed_bytes: u32, l: usize) -> String {
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
        tmp = tmp >> 2;
    }
    return res
}

// ---------------------
//      PrefixIter
// ---------------------

impl<'a> PrefixIter<'a> {
    pub fn new(stb: &SingleThreadedBuilder, packed_bytes: &'a [u128]) -> Self {
        let packed_u32_bytes = transform_to_u32_slice_unsafe(packed_bytes);
        Self {
            packed_u32_bytes: packed_u32_bytes,
            pos: 0,
            prefix_size: stb.prefix_size,
            current_u64_window: get_u64_unaligned(packed_u32_bytes, 0),
            remaining_inside_window: 16 - stb.prefix_size + 1,
            shift_inside_window: 0,
            filter_prefix_mask: std::u32::MAX >> (32 - 2*stb.prefix_size),
        }
    }
}

impl<'a> Iterator for PrefixIter<'a> {
    type Item = u32;
    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining_inside_window == 0 {
            if self.pos == self.packed_u32_bytes.len() - 2 {
                return None
            }

            self.pos += 1;
            self.current_u64_window = get_u64_unaligned(self.packed_u32_bytes, self.pos);
            self.remaining_inside_window = 16 - self.prefix_size + 1;
            self.shift_inside_window += 0;
        }

        let current = ((self.current_u64_window >> self.shift_inside_window) as u32) & self.filter_prefix_mask;
        self.shift_inside_window += 2;
        self.remaining_inside_window -= 1;
        return Some(current)
    }
}


// ---------------------
//      SuffixIter
// ---------------------

struct SuffixIter {

}

// --------------------------------
//      SingleThreadedBuilder
// --------------------------------

impl<'a> SingleThreadedBuilder {
    pub fn new(skb: &SketchBuilder, masks: Vec<u32>) -> Self {
        return Self {
            k: skb.k,
            prefix_size: skb.prefix_size,
            masks: masks,
        }
    }

    fn get_prefix_iterator(&self, packed_bytes: &'a [u128]) -> PrefixIter<'a> {
       return  PrefixIter::new(&self, packed_bytes)
    }

    pub fn build(&self, packed_bytes: &'a [u128]) {
        for kmer in self.get_prefix_iterator(packed_bytes) {
            let seq_as_string = packed_to_string(kmer, self.prefix_size);
            println!("{seq_as_string}");
        }
    }
}

// ------------------------
//      SketchBuilder
// ------------------------

impl SketchBuilder {
    pub fn new(k: usize, prefix_size: usize, nb_threads: usize) -> Self {
        Self {
            k,
            prefix_size,
            nb_threads,
        }
    }

    pub fn build(&self, seq: &PackedDNA) -> LexicSketch {
        let masks = create_masks(self);
        let (packed_bytes, _) = seq.bits();
        let kmer_prefix_mask: u32 = std::u32::MAX >> (32 - 2*self.prefix_size);
        let kmer_suffix_mask: u32 = std::u32::MAX >> (32 - 2*(self.k - self.prefix_size));
        let single_thread_builder = SingleThreadedBuilder::new(self, masks);
        single_thread_builder.build(packed_bytes);
        todo!()
    }
}
