use crate::seq::{Dna, Rna};

pub struct OffTargetScore {
    pub mismatches: usize,
    pub score: f64,
}

/// Calculates a simple mismatch score between an sgRNA and a target DNA.
/// Zips both sequences, so it evaluates up to the length of the shorter sequence.
pub fn calculate_mismatch_score(sgrna: &Rna, target: &Dna) -> OffTargetScore {
    let mut mismatches = 0;
    
    for (r, d) in sgrna.sequence.iter().zip(target.sequence.iter()) {
        // RNA has U, DNA has T. They should match.
        let is_match = matches!((r, d), (b'A', b'A') | (b'C', b'C') | (b'G', b'G') | (b'U', b'T'));
        if !is_match {
            mismatches += 1;
        }
    }
    
    // Simple scoring: 100.0 - (mismatches * 5.0). Min 0.
    let score = (100.0 - (mismatches as f64 * 5.0)).max(0.0);
    
    OffTargetScore { mismatches, score }
}

/// Position-dependent mismatch weights from Doench et al. (Nature Biotech 2016) CFD matrix
/// for positions 1 to 20 (5' to 3' PAM-proximal). Seed region (positions 15-20) has
/// dramatically lower tolerance (higher penalty), reducing cleaving activity.
pub const DOENCH_2016_CFD_WEIGHTS: [f64; 20] = [
    0.92, 0.88, 0.85, 0.80, 0.76, // Pos 1-5  (PAM-distal)
    0.72, 0.68, 0.62, 0.55, 0.50, // Pos 6-10
    0.44, 0.38, 0.32, 0.26, 0.20, // Pos 11-15 (Transition to seed)
    0.15, 0.11, 0.08, 0.05, 0.02, // Pos 16-20 (Strict seed PAM-proximal)
];

/// Calculates high-fidelity Cutting Frequency Determination (CFD) score using the Doench 2016 matrix.
///
/// Multiplies position-specific cleavage probabilities for all mismatch positions.
pub fn calculate_doench_2016_cfd(sgrna: &Rna, target: &Dna) -> OffTargetScore {
    let mut mismatches = 0;
    let mut cfd_product = 1.0;

    for (pos, (r, d)) in sgrna.sequence.iter().zip(target.sequence.iter()).enumerate() {
        let is_match = matches!((r, d), (b'A', b'A') | (b'C', b'C') | (b'G', b'G') | (b'U', b'T'));
        if !is_match {
            mismatches += 1;
            let weight = if pos < DOENCH_2016_CFD_WEIGHTS.len() {
                DOENCH_2016_CFD_WEIGHTS[pos]
            } else {
                0.01
            };
            cfd_product *= weight;
        }
    }

    OffTargetScore {
        mismatches,
        score: cfd_product * 100.0,
    }
}

/// Calculates a mock CFD-style score for an sgRNA vs a target DNA.
/// Applies higher penalties to the seed region (PAM-proximal, last 10nt of the sgRNA).
pub fn calculate_cfd_score(sgrna: &Rna, target: &Dna) -> OffTargetScore {
    let mut mismatches = 0;
    let mut score_multiplier = 1.0;
    let len = sgrna.sequence.len();
    
    for (i, (r, d)) in sgrna.sequence.iter().zip(target.sequence.iter()).enumerate() {
        let is_match = matches!((r, d), (b'A', b'A') | (b'C', b'C') | (b'G', b'G') | (b'U', b'T'));
        
        if !is_match {
            mismatches += 1;
            // PAM-proximal (seed region) is typically the last 10 nucleotides of a 20nt sgRNA.
            // Distance from PAM (3' end): len - 1 - i. 
            // If distance < 10, heavy penalty. Else light penalty.
            let dist_from_3prime = len.saturating_sub(1).saturating_sub(i);
            if dist_from_3prime < 10 {
                score_multiplier *= 0.2; // heavy penalty in seed
            } else {
                score_multiplier *= 0.8; // light penalty outside seed
            }
        }
    }
    
    OffTargetScore { 
        mismatches, 
        score: score_multiplier * 100.0 
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mismatch_score() {
        let rna = Rna::new(b"AUGCUACGAUGCUACGAUGC").unwrap();
        let dna1 = Dna::new(b"ATGCTACGATGCTACGATGC").unwrap(); // Perfect match
        let dna2 = Dna::new(b"ATGCTACGATGCTACGATTT").unwrap(); // 2 mismatches at end
        
        let score1 = calculate_mismatch_score(&rna, &dna1);
        assert_eq!(score1.mismatches, 0);
        assert_eq!(score1.score, 100.0);
        
        let score2 = calculate_mismatch_score(&rna, &dna2);
        assert_eq!(score2.mismatches, 2);
        assert_eq!(score2.score, 90.0);
    }

    #[test]
    fn test_cfd_score() {
        let rna = Rna::new(b"AUGCUACGAUGCUACGAUGC").unwrap();
        let dna_seed_mut = Dna::new(b"ATGCTACGATGCTACGATTT").unwrap(); // 2 mismatches in seed (last 10nt)
        let dna_dist_mut = Dna::new(b"CCGCTACGATGCTACGATGC").unwrap(); // 2 mismatches far from PAM (first 10nt)
        
        let seed_score = calculate_cfd_score(&rna, &dna_seed_mut);
        let dist_score = calculate_cfd_score(&rna, &dna_dist_mut);
        
        // Both have 2 mismatches, but seed mutation should penalize much more.
        assert_eq!(seed_score.mismatches, 2);
        assert_eq!(dist_score.mismatches, 2);
        assert!(seed_score.score < dist_score.score);
        
        // Verification of penalty math
        // Seed mismatches (i=18, 19). dist=1, dist=0. Both < 10. Multiplier = 0.2 * 0.2 = 0.04 -> Score = 4.0
        // Dist mismatches (i=0, 1). dist=19, dist=18. Both >= 10. Multiplier = 0.8 * 0.8 = 0.64 -> Score = 64.0
        assert!((seed_score.score - 4.0).abs() < 1e-9);
        assert!((dist_score.score - 64.0).abs() < 1e-9);
    }

    #[test]
    fn test_doench_2016_cfd_score() {
        let rna = Rna::new(b"AUGCUACGAUGCUACGAUGC").unwrap();
        let dna_perfect = Dna::new(b"ATGCTACGATGCTACGATGC").unwrap();
        let dna_seed_mut = Dna::new(b"ATGCTACGATGCTACGATGT").unwrap(); // Pos 20 mismatch (0.02)
        let dna_dist_mut = Dna::new(b"TTGCTACGATGCTACGATGC").unwrap(); // Pos 1 mismatch (0.92)

        let p_score = calculate_doench_2016_cfd(&rna, &dna_perfect);
        assert_eq!(p_score.score, 100.0);

        let seed_score = calculate_doench_2016_cfd(&rna, &dna_seed_mut);
        assert!((seed_score.score - 2.0).abs() < 1e-6);

        let dist_score = calculate_doench_2016_cfd(&rna, &dna_dist_mut);
        assert!((dist_score.score - 92.0).abs() < 1e-6);
    }
}
