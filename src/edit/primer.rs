use crate::seq::Dna;
use crate::align::SmithWaterman;

pub struct Primer {
    pub sequence: Dna,
    pub tm: f64,
    pub gc_percent: f64,
}

impl Primer {
    pub fn new(sequence: Dna) -> Self {
        let tm = calculate_tm(&sequence);
        let gc_percent = sequence.gc_content();
        Self { sequence, tm, gc_percent }
    }
}

/// Calculates the melting temperature (Tm) of a DNA sequence.
/// Uses the Wallace rule for sequences < 14nt (2*(A+T) + 4*(G+C))
/// Uses a simplified salt-adjusted formula for sequences >= 14nt:
/// Tm = 64.9 + 41.0 * (num_G + num_C - 16.4) / length
pub fn calculate_tm(primer: &Dna) -> f64 {
    let bytes = &primer.sequence;
    let len = bytes.len() as f64;
    if len == 0.0 {
        return 0.0;
    }
    
    let mut num_g = 0;
    let mut num_c = 0;
    let mut num_a = 0;
    let mut num_t = 0;
    
    for &b in bytes {
        match b {
            b'G' => num_g += 1,
            b'C' => num_c += 1,
            b'A' => num_a += 1,
            b'T' => num_t += 1,
            _ => {}
        }
    }
    
    if len < 14.0 {
        (2 * (num_a + num_t) + 4 * (num_g + num_c)) as f64
    } else {
        64.9 + 41.0 * ((num_g + num_c) as f64 - 16.4) / len
    }
}

/// Checks if a primer has a GC clamp.
/// A simple definition used here is the presence of G or C at the 3' end (the last nucleotide).
pub fn check_gc_clamp(primer: &Dna) -> bool {
    if let Some(&last) = primer.sequence.last() {
        last == b'G' || last == b'C'
    } else {
        false
    }
}

/// Detects if a primer is likely to form a hairpin.
/// A simplified check: performs local alignment of the primer against its own reverse complement.
/// High score implies a stable stem-loop formation.
pub fn detect_hairpin(primer: &Dna) -> bool {
    let rc = primer.reverse_complement();
    let sw = SmithWaterman::new(2, -2, -2);
    let res = sw.align(&primer.sequence, &rc.sequence);
    
    // A score >= 8 means at least 4 perfectly matching bases.
    res.score >= 8
}

/// Detects if two primers are likely to form a primer dimer.
/// Performs local alignment of primer1 against the reverse complement of primer2.
pub fn detect_dimer(primer1: &Dna, primer2: &Dna) -> bool {
    let rc2 = primer2.reverse_complement();
    let sw = SmithWaterman::new(2, -2, -2);
    let res = sw.align(&primer1.sequence, &rc2.sequence);
    
    // A score >= 10 means at least 5 perfectly matching bases.
    res.score >= 10
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_tm() {
        let short_dna = Dna::new(b"ATGC").unwrap(); // 2*(2) + 4*(2) = 12
        assert_eq!(calculate_tm(&short_dna), 12.0);
        
        let long_dna = Dna::new(b"ATGCATGCATGCATGC").unwrap(); // length 16, GC 8
        let expected_tm = 64.9 + 41.0 * (8.0 - 16.4) / 16.0;
        assert!((calculate_tm(&long_dna) - expected_tm).abs() < 1e-9);
    }

    #[test]
    fn test_gc_clamp() {
        let clamp_dna = Dna::new(b"ATGC").unwrap(); // Ends in C
        assert!(check_gc_clamp(&clamp_dna));
        
        let no_clamp_dna = Dna::new(b"ATGA").unwrap(); // Ends in A
        assert!(!check_gc_clamp(&no_clamp_dna));
    }

    #[test]
    fn test_hairpin() {
        let hairpin_dna = Dna::new(b"GCATGC").unwrap(); 
        assert!(detect_hairpin(&hairpin_dna));
        
        let no_hairpin_dna = Dna::new(b"AAAAAAAA").unwrap(); // RC = TTTTTTTT (score 0)
        assert!(!detect_hairpin(&no_hairpin_dna));
    }

    #[test]
    fn test_dimer() {
        let p1 = Dna::new(b"AAAAACCC").unwrap();
        let p2 = Dna::new(b"GGGTTTTT").unwrap(); // RC of p2 is AAAAACCC
        assert!(detect_dimer(&p1, &p2));
        
        let p3 = Dna::new(b"AAAAACCC").unwrap();
        let p4 = Dna::new(b"AAAAACCC").unwrap(); // RC is GGGTTTTT (score 0)
        assert!(!detect_dimer(&p3, &p4));
    }
}
