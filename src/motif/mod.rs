use crate::seq::Dna;

/// A Position Weight Matrix (PWM) for scanning motifs in DNA sequences.
/// Columns represent positions; rows represent A, C, G, T respectively.
#[derive(Debug, Clone)]
pub struct Pwm {
    pub weights: Vec<[f64; 4]>, // weights[position][ACGT]
    pub threshold: f64,
}

fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

impl Pwm {
    /// Creates a PWM from a counts matrix (rows=positions, cols=ACGT) with pseudocount.
    pub fn from_counts(counts: &[[f64; 4]], pseudocount: f64) -> Self {
        let weights = counts.iter().map(|row| {
            let total: f64 = row.iter().sum::<f64>() + 4.0 * pseudocount;
            let freqs = [
                (row[0] + pseudocount) / total,
                (row[1] + pseudocount) / total,
                (row[2] + pseudocount) / total,
                (row[3] + pseudocount) / total,
            ];
            // Log-odds against background (uniform 0.25)
            [
                freqs[0].ln() - 0.25_f64.ln(),
                freqs[1].ln() - 0.25_f64.ln(),
                freqs[2].ln() - 0.25_f64.ln(),
                freqs[3].ln() - 0.25_f64.ln(),
            ]
        }).collect();

        Self { weights, threshold: 0.0 }
    }

    /// Returns the number of positions in this PWM.
    pub fn len(&self) -> usize {
        self.weights.len()
    }

    /// Scores a DNA k-mer against this PWM. Returns None for invalid characters.
    pub fn score(&self, kmer: &[u8]) -> Option<f64> {
        if kmer.len() != self.weights.len() {
            return None;
        }
        let mut total = 0.0;
        for (i, &b) in kmer.iter().enumerate() {
            let idx = base_index(b)?;
            total += self.weights[i][idx];
        }
        Some(total)
    }

    /// Scans a DNA sequence and returns (position, score) for all windows scoring >= threshold.
    pub fn scan(&self, seq: &Dna) -> Vec<(usize, f64)> {
        let motif_len = self.weights.len();
        let seq_len = seq.sequence.len();
        if motif_len > seq_len {
            return Vec::new();
        }
        let mut hits = Vec::new();
        for i in 0..=(seq_len - motif_len) {
            if let Some(s) = self.score(&seq.sequence[i..i + motif_len]) {
                if s >= self.threshold {
                    hits.push((i, s));
                }
            }
        }
        hits
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_pwm() -> Pwm {
        // A simple motif favouring "ACGT"
        let counts = [
            [10.0, 0.0, 0.0, 0.0],  // pos 0 -> A
            [0.0, 10.0, 0.0, 0.0],  // pos 1 -> C
            [0.0, 0.0, 10.0, 0.0],  // pos 2 -> G
            [0.0, 0.0, 0.0, 10.0],  // pos 3 -> T
        ];
        let mut pwm = Pwm::from_counts(&counts, 0.1);
        pwm.threshold = 5.0;
        pwm
    }

    #[test]
    fn test_pwm_score_perfect_match() {
        let pwm = make_pwm();
        let score = pwm.score(b"ACGT").unwrap();
        assert!(score > 5.0);
    }

    #[test]
    fn test_pwm_score_mismatch() {
        let pwm = make_pwm();
        let score = pwm.score(b"TTTT").unwrap();
        assert!(score < 0.0);
    }

    #[test]
    fn test_pwm_scan() {
        let pwm = make_pwm();
        let dna = Dna::new(b"NNNNACGTNNNN").unwrap();
        let hits = pwm.scan(&dna);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].0, 4); // ACGT starts at index 4
    }
}
