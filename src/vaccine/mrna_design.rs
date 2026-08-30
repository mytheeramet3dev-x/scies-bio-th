//! Full-length mRNA cancer vaccine construct assembly and human codon optimization.
//!
//! Assembles a clinical-grade synthetic mRNA sequence:
//! `[5' Cap] -> [5' UTR] -> [Kozak] -> [Signal Peptide] -> [Neoepitope Polytope (with AAY linkers)] -> [MITD] -> [Dual Stop] -> [3' UTR] -> [Poly-A Tail]`

use crate::error::{BioError, Result};
use crate::seq::{Dna, Protein, Rna};
use crate::vaccine::neoantigen::Neoantigen;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Human alpha-globin optimized 5' UTR (35 nt).
pub const HUMAN_5_PRIME_UTR: &[u8] = b"TCTTCTGGTCCCCACAGACTCAGAGAGAACCCACC";

/// Kozak consensus sequence for maximal translation initiation.
pub const KOZAK_CONSENSUS: &[u8] = b"GCCACC";

/// Sec secretion signal peptide (20 aa) directing nascent polypeptide to the ER.
pub const SECRETION_SIGNAL_PEPTIDE: &[u8] = b"MYRMQLLSCIALSLALVTNS";

/// Proteasome cleavage linker peptide (`AAY`) promoting efficient epitope excision in dendritic cells.
pub const PROTEASOME_LINKER_AAY: &[u8] = b"AAY";

/// MHC Class I Trafficking Domain (MITD, 35 aa) directing antigens to endosomal/lysosomal compartments.
pub const MITD_TRAFFICKING_DOMAIN: &[u8] = b"IVGIVAGLAVLAVVVIGAVVATVMCRRKSSGGKGGSYSQAAS";

/// Dual human beta-globin 3' UTR enhancing mRNA cytoplasmic stability.
pub const DUAL_3_PRIME_UTR: &[u8] =
    b"GCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGCA";

/// Assembled clinical-grade personalized mRNA vaccine construct.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MrnaVaccineConstruct {
    /// Patient or cohort identifier.
    pub patient_id: String,
    /// List of targeted gene neoantigens included in the polytope.
    pub included_genes: Vec<String>,
    /// Full open reading frame (ORF) translated protein sequence.
    pub polytope_protein: Protein,
    /// Complete synthetic mRNA sequence (including UTRs and 120-nt Poly-A tail).
    pub full_mrna_transcript: Rna,
    /// Total transcript length in nucleotides.
    pub total_length_nt: usize,
    /// GC content percentage of the codon-optimized ORF.
    pub orf_gc_content: f64,
}

impl MrnaVaccineConstruct {
    /// Assembles a personalized multi-epitope mRNA vaccine from a prioritized list of neoantigens.
    pub fn build(patient_id: &str, neoantigens: &[Neoantigen]) -> Result<Self> {
        if neoantigens.is_empty() {
            return Err(BioError::ParseError(
                "Cannot build mRNA vaccine without neoantigens".to_string(),
            ));
        }

        let mut included_genes = Vec::new();
        let mut orf_amino_acids = Vec::new();

        // 1. Prepend Sec Signal Peptide
        orf_amino_acids.extend_from_slice(SECRETION_SIGNAL_PEPTIDE);

        // 2. Concatenate Neoantigen Peptides with AAY Proteasome Cleavage Linkers
        for (i, neo) in neoantigens.iter().enumerate() {
            if i > 0 {
                orf_amino_acids.extend_from_slice(PROTEASOME_LINKER_AAY);
            }
            orf_amino_acids.extend_from_slice(&neo.mut_peptide.sequence);
            included_genes.push(format!("{}_{}{}{}", neo.gene_symbol, neo.wt_aa, neo.original_protein_pos, neo.mut_aa));
        }

        // 3. Append MITD Trafficking Domain
        orf_amino_acids.extend_from_slice(PROTEASOME_LINKER_AAY);
        orf_amino_acids.extend_from_slice(MITD_TRAFFICKING_DOMAIN);

        let polytope_protein = Protein::new(&orf_amino_acids)?;

        // 4. Perform Human Codon Optimization on the ORF
        let optimized_orf_dna = human_codon_optimize(&polytope_protein)?;
        let orf_gc = optimized_orf_dna.gc_content();

        // 5. Assemble Full mRNA Sequence: 5' UTR + Kozak + ORF + Stop + 3' UTR + PolyA(120)
        let mut full_dna_builder = Vec::new();
        full_dna_builder.extend_from_slice(HUMAN_5_PRIME_UTR);
        full_dna_builder.extend_from_slice(KOZAK_CONSENSUS);
        full_dna_builder.extend_from_slice(&optimized_orf_dna.sequence);
        full_dna_builder.extend_from_slice(b"TGATAA"); // Tandem dual stop codons
        full_dna_builder.extend_from_slice(DUAL_3_PRIME_UTR);
        full_dna_builder.extend(std::iter::repeat(b'A').take(120)); // 120-nt Poly-A tail

        // Convert DNA representation to mature mRNA (T -> U)
        let mut rna_bytes = Vec::with_capacity(full_dna_builder.len());
        for &b in &full_dna_builder {
            if b == b'T' {
                rna_bytes.push(b'U');
            } else {
                rna_bytes.push(b);
            }
        }
        let full_mrna_transcript = Rna::new(&rna_bytes)?;
        let total_length_nt = full_mrna_transcript.sequence.len();

        Ok(Self {
            patient_id: patient_id.to_string(),
            included_genes,
            polytope_protein,
            full_mrna_transcript,
            total_length_nt,
            orf_gc_content: orf_gc,
        })
    }
}

/// Optimizes an amino acid sequence using highest-frequency human codons and GC enrichment.
pub fn human_codon_optimize(protein: &Protein) -> Result<Dna> {
    let mut dna_bytes = Vec::with_capacity(protein.sequence.len() * 3);

    for &aa in &protein.sequence {
        let codon: &[u8] = match aa {
            b'A' => b"GCC", // Alanine (40% human usage, high GC)
            b'C' => b"TGC", // Cysteine (55%)
            b'D' => b"GAC", // Aspartic Acid (54%)
            b'E' => b"GAG", // Glutamic Acid (58%)
            b'F' => b"TTC", // Phenylalanine (58%)
            b'G' => b"GGC", // Glycine (34%)
            b'H' => b"CAC", // Histidine (59%)
            b'I' => b"ATC", // Isoleucine (48%)
            b'K' => b"AAG", // Lysine (58%)
            b'L' => b"CTG", // Leucine (41%)
            b'M' => b"ATG", // Methionine (100%)
            b'N' => b"AAC", // Asparagine (53%)
            b'P' => b"CCC", // Proline (33%)
            b'Q' => b"CAG", // Glutamine (75%)
            b'R' => b"CGG", // Arginine (21%)
            b'S' => b"AGC", // Serine (24%)
            b'T' => b"ACC", // Threonine (36%)
            b'V' => b"GTG", // Valine (47%)
            b'W' => b"TGG", // Tryptophan (100%)
            b'Y' => b"TAC", // Tyrosine (57%)
            _ => b"GCC",   // Default fallback
        };
        dna_bytes.extend_from_slice(codon);
    }

    Dna::new(&dna_bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mrna_vaccine_construct_assembly() {
        let kras_wt = Protein::new(b"MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSY").unwrap();
        let kras_neo = Neoantigen::from_mutation("KRAS", &kras_wt, 12, 'D', 6).unwrap();

        let tp53_wt = Protein::new(b"MSSSRVPQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNK").unwrap();
        let tp53_neo = Neoantigen::from_mutation("TP53", &tp53_wt, 15, 'H', 6).unwrap();

        let vaccine = MrnaVaccineConstruct::build("PATIENT_001_MELANOMA", &[kras_neo, tp53_neo]).unwrap();

        assert_eq!(vaccine.patient_id, "PATIENT_001_MELANOMA");
        assert_eq!(vaccine.included_genes.len(), 2);
        assert!(vaccine.orf_gc_content > 55.0); // High GC from human codon optimization
        assert!(vaccine.total_length_nt > 300);

        // Verify Poly-A tail at 3' end
        let end_slice = &vaccine.full_mrna_transcript.sequence[vaccine.total_length_nt - 50..];
        for &b in end_slice {
            assert_eq!(b, b'A');
        }
    }
}
