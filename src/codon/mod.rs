use crate::seq::{Dna, Protein, Rna};
use std::collections::HashMap;

pub struct CodonTable {
    table: HashMap<[u8; 3], u8>,
}

impl CodonTable {
    /// Creates the standard genetic code codon table.
    pub fn standard() -> Self {
        let mut table = HashMap::new();
        let mappings = [
            (b"UUU", b'F'), (b"UUC", b'F'), (b"UUA", b'L'), (b"UUG", b'L'),
            (b"CUU", b'L'), (b"CUC", b'L'), (b"CUA", b'L'), (b"CUG", b'L'),
            (b"AUU", b'I'), (b"AUC", b'I'), (b"AUA", b'I'), (b"AUG", b'M'),
            (b"GUU", b'V'), (b"GUC", b'V'), (b"GUA", b'V'), (b"GUG", b'V'),
            
            (b"UCU", b'S'), (b"UCC", b'S'), (b"UCA", b'S'), (b"UCG", b'S'),
            (b"CCU", b'P'), (b"CCC", b'P'), (b"CCA", b'P'), (b"CCG", b'P'),
            (b"ACU", b'T'), (b"ACC", b'T'), (b"ACA", b'T'), (b"ACG", b'T'),
            (b"GCU", b'A'), (b"GCC", b'A'), (b"GCA", b'A'), (b"GCG", b'A'),
            
            (b"UAU", b'Y'), (b"UAC", b'Y'), (b"UAA", b'*'), (b"UAG", b'*'),
            (b"CAU", b'H'), (b"CAC", b'H'), (b"CAA", b'Q'), (b"CAG", b'Q'),
            (b"AAU", b'N'), (b"AAC", b'N'), (b"AAA", b'K'), (b"AAG", b'K'),
            (b"GAU", b'D'), (b"GAC", b'D'), (b"GAA", b'E'), (b"GAG", b'E'),
            
            (b"UGU", b'C'), (b"UGC", b'C'), (b"UGA", b'*'), (b"UGG", b'W'),
            (b"CGU", b'R'), (b"CGC", b'R'), (b"CGA", b'R'), (b"CGG", b'R'),
            (b"AGU", b'S'), (b"AGC", b'S'), (b"AGA", b'R'), (b"AGG", b'R'),
            (b"GGU", b'G'), (b"GGC", b'G'), (b"GGA", b'G'), (b"GGG", b'G'),
        ];
        
        for (codon, aa) in mappings {
            let mut arr = [0u8; 3];
            arr.copy_from_slice(codon);
            table.insert(arr, aa);
        }
        
        Self { table }
    }

    /// Translates an RNA sequence into a Protein sequence using this codon table.
    pub fn translate(&self, rna: &Rna) -> crate::Result<Protein> {
        let mut protein_seq = Vec::new();
        let bytes = &rna.sequence;
        
        for chunk in bytes.chunks_exact(3) {
            let mut codon = [0u8; 3];
            codon.copy_from_slice(chunk);
            if let Some(&aa) = self.table.get(&codon) {
                protein_seq.push(aa);
            } else {
                protein_seq.push(b'X'); // Unknown amino acid
            }
        }
        
        Protein::new(&protein_seq)
    }
}

pub struct OrfFinder {
    pub min_len_codons: usize,
}

impl OrfFinder {
    /// Creates a new OrfFinder with a default minimum length of 30 codons (90 bp).
    pub fn new() -> Self {
        Self { min_len_codons: 30 }
    }

    /// Finds Open Reading Frames (ORFs) in a DNA sequence.
    /// An ORF starts with ATG and ends with a stop codon (TAA, TAG, TGA).
    /// Scans all 6 frames (3 forward, 3 reverse complement).
    pub fn find_orfs(&self, seq: &Dna) -> Vec<Dna> {
        let mut orfs = Vec::new();
        
        let frames = [
            seq.clone(),
            seq.reverse_complement(),
        ];
        
        for frame in frames {
            let bytes = &frame.sequence;
            for offset in 0..3 {
                let mut start_idx = None;
                let mut i = offset;
                while i + 2 < bytes.len() {
                    let codon = &bytes[i..i+3];
                    if codon == b"ATG" && start_idx.is_none() {
                        start_idx = Some(i);
                    } else if codon == b"TAA" || codon == b"TAG" || codon == b"TGA" {
                        if let Some(start) = start_idx {
                            let len_codons = (i + 3 - start) / 3;
                            if len_codons >= self.min_len_codons {
                                if let Ok(orf_dna) = Dna::new(&bytes[start..i+3]) {
                                    orfs.push(orf_dna);
                                }
                            }
                            start_idx = None;
                        }
                    }
                    i += 3;
                }
            }
        }
        
        orfs
    }
}

impl Default for OrfFinder {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_codon_translation() {
        let table = CodonTable::standard();
        let rna = Rna::new(b"AUGGCCUGGUAA").unwrap(); // M A W *
        let protein = table.translate(&rna).unwrap();
        assert_eq!(protein.sequence, b"MAW*");
    }

    #[test]
    fn test_orf_finder() {
        let mut finder = OrfFinder::new();
        finder.min_len_codons = 3; // min 3 codons (9 bp) including stop
        
        // Forward frame 0: ATG CCC TAA (9bp) -> valid
        let dna = Dna::new(b"ATGCCCTAAGGG").unwrap();
        let orfs = finder.find_orfs(&dna);
        assert_eq!(orfs.len(), 1);
        assert_eq!(orfs[0].sequence, b"ATGCCCTAA");
    }

    #[test]
    fn test_orf_finder_reverse() {
        let mut finder = OrfFinder::new();
        finder.min_len_codons = 3;
        
        // Reverse complement of ATG CCC TAA is TTA GGG CAT
        // If we put TTAGGGCAT, its reverse complement is ATGCCCTAA which has a valid ORF
        let dna = Dna::new(b"TTAGGGCAT").unwrap();
        let orfs = finder.find_orfs(&dna);
        assert_eq!(orfs.len(), 1);
        assert_eq!(orfs[0].sequence, b"ATGCCCTAA");
    }
}
