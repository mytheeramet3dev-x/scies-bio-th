use crate::seq::{Dna, Rna};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PamSite {
    pub index: usize,
    pub sequence: Vec<u8>,
}

/// Finds SpCas9 PAM sites (NGG) in a DNA sequence.
pub fn find_spcas9_pam(seq: &Dna) -> Vec<PamSite> {
    let mut sites = Vec::new();
    let bytes = &seq.sequence;
    if bytes.len() < 3 {
        return sites;
    }
    
    for i in 0..=bytes.len() - 3 {
        let n = bytes[i];
        let g1 = bytes[i + 1];
        let g2 = bytes[i + 2];
        
        if (n == b'A' || n == b'C' || n == b'G' || n == b'T' || n == b'N') && g1 == b'G' && g2 == b'G' {
            sites.push(PamSite {
                index: i,
                sequence: vec![n, g1, g2],
            });
        }
    }
    sites
}

/// Finds Cas12a PAM sites (TTTV, where V is A, C, or G) in a DNA sequence.
pub fn find_cas12a_pam(seq: &Dna) -> Vec<PamSite> {
    let mut sites = Vec::new();
    let bytes = &seq.sequence;
    if bytes.len() < 4 {
        return sites;
    }
    
    for i in 0..=bytes.len() - 4 {
        let t1 = bytes[i];
        let t2 = bytes[i + 1];
        let t3 = bytes[i + 2];
        let v = bytes[i + 3];
        
        if t1 == b'T' && t2 == b'T' && t3 == b'T' && (v == b'A' || v == b'C' || v == b'G') {
            sites.push(PamSite {
                index: i,
                sequence: vec![t1, t2, t3, v],
            });
        }
    }
    sites
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SgRnaCandidate {
    pub sequence: Rna,
    pub pam_index: usize,
}

/// Generates sgRNA candidates for SpCas9 (20nt before the NGG PAM).
pub fn generate_sgrna_candidates(seq: &Dna) -> Vec<SgRnaCandidate> {
    let mut candidates = Vec::new();
    let pams = find_spcas9_pam(seq);
    
    for pam in pams {
        if pam.index >= 20 {
            let start = pam.index - 20;
            let dna_slice = &seq.sequence[start..pam.index];
            // Convert DNA slice to RNA (T -> U)
            let mut rna_seq = Vec::with_capacity(20);
            for &b in dna_slice {
                if b == b'T' {
                    rna_seq.push(b'U');
                } else {
                    rna_seq.push(b);
                }
            }
            if let Ok(rna) = Rna::new(&rna_seq) {
                candidates.push(SgRnaCandidate {
                    sequence: rna,
                    pam_index: pam.index,
                });
            }
        }
    }
    candidates
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_spcas9_pam() {
        let dna = Dna::new(b"ATGCAGGTACTGG").unwrap();
        let pams = find_spcas9_pam(&dna);
        assert_eq!(pams.len(), 2);
        assert_eq!(pams[0].index, 4); // AGG
        assert_eq!(pams[1].index, 10); // TGG
    }

    #[test]
    fn test_find_cas12a_pam() {
        let dna = Dna::new(b"TTTACGTTTGG").unwrap();
        let pams = find_cas12a_pam(&dna);
        assert_eq!(pams.len(), 2);
        assert_eq!(pams[0].index, 0); // TTTA
        assert_eq!(pams[1].index, 6); // TTTG
    }

    #[test]
    fn test_sgrna_candidates() {
        let seq = b"ATGCTACGATGCTACGATGCAGGT"; // 24 bp, PAM at 20 (AGGT -> NGG is AGG at 20)
        let dna = Dna::new(seq).unwrap();
        let candidates = generate_sgrna_candidates(&dna);
        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].pam_index, 20);
        // RNA equivalent of "ATGCTACGATGCTACGATGC" (first 20 chars) -> "AUGCUACGAUGCUACGAUGC"
        assert_eq!(candidates[0].sequence.sequence, b"AUGCUACGAUGCUACGAUGC");
    }
}
