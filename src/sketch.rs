pub struct LexicSketch {
    k: u8,
    prefix_size: u8,
    fingerprint: Vec<u32>
}

impl LexicSketch {
    pub fn get_k(&self) -> u8 {
        self.k
    }

    pub fn get_prefix_size(&self) -> u8 {
        self.prefix_size
    }

    pub fn get_fingerprint(&self) -> &[u32] {
        &self.fingerprint
    }

    pub fn deserialize() {
        // TODO
    }

    pub fn compare(&self, sk: &Self) -> Vec<u32> {
        // TODO
        return vec![]
    }

    pub fn get_score(&self, sk: &Self) -> f64 {
        // TODO
        return 0.
    }
}
