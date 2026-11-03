//! DSB-free precision genome editing: Base Editing (ABE8e / BE4max), Prime Editing (pegRNA),
//! and Multiplex Chromosomal Translocation risk modeling.
//!
//! References:
//! - Gaudelli et al. (2017) *Nature* 551:464 / Richter et al. (2020) *Nat Biotechnol* 38:884 (ABE8e).
//! - Komor et al. (2016) *Nature* 533:420 / Koblan et al. (2018) *Nat Biotechnol* 36:843 (BE4max).
//! - Anzalone et al. (2019) *Nature* 576:149 (Prime Editing pegRNA: Spacer + PBS + RTT).
//! - Bothmer et al. (2017) *Nat Commun* 8:13905 (Multiplex DSB chromosomal translocation kinetics).

use crate::edit::crispr::find_spcas9_pam;
use crate::error::{BioError, Result};
use crate::seq::{Dna, Rna};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Class of DSB-free deaminase Base Editor fused to Cas9 D10A nickase (`nCas9`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum BaseEditorType {
    /// Adenine Base Editor (`ABE8e`, evolved TadA-8e): catalyzes $A \cdot T \rightarrow G \cdot C$
    /// within protospacer positions $4\text{--}8$ (1-indexed from PAM-distal 5' end).
    Abe8e,
    /// Cytosine Base Editor (`BE4max`, rAPOBEC1 + $2\times$ UGI): catalyzes $C \cdot G \rightarrow T \cdot A$
    /// within protospacer positions $4\text{--}8$ (1-indexed from PAM-distal 5' end).
    Be4Max,
}

impl BaseEditorType {
    /// Returns `(source_base, converted_base)` on the non-target protospacer DNA strand.
    pub fn transition_pair(&self) -> (u8, u8) {
        match self {
            BaseEditorType::Abe8e => (b'A', b'G'),
            BaseEditorType::Be4Max => (b'C', b'T'),
        }
    }

    /// Canonical high-efficiency deamination window `[start_pos_1idx, end_pos_1idx]`
    /// along the 20-nt protospacer (1 = PAM-distal 5' end, 20 = PAM-proximal 3' end).
    pub fn activity_window_1idx(&self) -> (usize, usize) {
        (4, 8)
    }
}

/// A candidate Base Editing guide placing a target nucleotide inside the deamination window.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct BaseEditCandidate {
    /// Base editor architecture (`Abe8e` or `Be4Max`).
    pub editor: BaseEditorType,
    /// 20-nt sgRNA protospacer sequence.
    pub spacer: Rna,
    /// 0-based index of the `NGG` PAM in the input DNA sequence.
    pub pam_index: usize,
    /// 1-based position ($4\text{--}8$) of the target base within the 20-nt protospacer.
    pub target_protospacer_pos_1idx: usize,
    /// Number of unintended bystander bases of the same nucleotide type inside the $4\text{--}8$ window.
    pub bystander_count: usize,
    /// Predicted product purity and editing efficiency score ($0.0$ to $1.0$), penalized by bystanders.
    pub precision_score: f64,
}

/// Scans `locus_dna` for SpCas9 (`NGG`) guides that position `target_index` (0-based coordinate in `locus_dna`)
/// within the canonical $4\text{--}8$ activity window of `editor`, counting bystander substrates.
pub fn design_base_editor_guides(
    locus_dna: &Dna,
    target_index: usize,
    editor: BaseEditorType,
) -> Result<Vec<BaseEditCandidate>> {
    let bytes = &locus_dna.sequence;
    if target_index >= bytes.len() {
        return Err(BioError::InvalidInput(format!(
            "target_index ({target_index}) is out of bounds for locus length ({})",
            bytes.len()
        )));
    }

    let (src_base, _) = editor.transition_pair();
    if bytes[target_index] != src_base {
        return Err(BioError::InvalidInput(format!(
            "Target base '{}' at index {} does not match {:?} substrate '{}'",
            bytes[target_index] as char, target_index, editor, src_base as char
        )));
    }

    let (win_start, win_end) = editor.activity_window_1idx();
    let pams = find_spcas9_pam(locus_dna);
    let mut candidates = Vec::new();

    for pam in pams {
        if pam.index < 20 {
            continue;
        }
        let proto_start = pam.index - 20;
        if target_index < proto_start || target_index >= pam.index {
            continue;
        }
        let pos_1idx = (target_index - proto_start) + 1;
        if pos_1idx < win_start || pos_1idx > win_end {
            continue;
        }

        let proto_slice = &bytes[proto_start..pam.index];
        // Count bystander substrates inside positions [win_start..=win_end] excluding target_index
        let mut bystanders = 0usize;
        for p in win_start..=win_end {
            let seq_idx = proto_start + (p - 1);
            if seq_idx != target_index && bytes[seq_idx] == src_base {
                bystanders += 1;
            }
        }

        // Optimal peak efficiency at position 6; penalized by bystander edits
        let center_dist = (pos_1idx as f64 - 6.0).abs();
        let positional_eff = (0.92 - center_dist * 0.08).clamp(0.50, 0.95);
        let bystander_penalty = 0.60_f64.powi(bystanders as i32);
        let precision_score = (positional_eff * bystander_penalty).clamp(0.05, 0.98);

        let rna_bytes: Vec<u8> = proto_slice
            .iter()
            .map(|&b| if b == b'T' { b'U' } else { b })
            .collect();
        if let Ok(spacer) = Rna::new(&rna_bytes) {
            candidates.push(BaseEditCandidate {
                editor,
                spacer,
                pam_index: pam.index,
                target_protospacer_pos_1idx: pos_1idx,
                bystander_count: bystanders,
                precision_score,
            });
        }
    }

    candidates.sort_by(|a, b| {
        b.precision_score
            .partial_cmp(&a.precision_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(candidates)
}

/// Prime Editing Guide RNA (`pegRNA`) construct for SpCas9(H840A) Reverse Transcriptase (`PE2`/`PE3`).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PegRnaDesign {
    /// 20-nt SpCas9 spacer sequence (5' -> 3').
    pub spacer: Rna,
    /// Primer Binding Site (`PBS`, 5' -> 3' RNA), reverse-complementary to the nicked PAM strand
    /// immediately 5' of the H840A single-strand nick site (`pam_index - 3`).
    pub pbs: Rna,
    /// Reverse Transcriptase Template (`RTT`, 5' -> 3' RNA), reverse-complementary to the desired
    /// edited DNA flap downstream of the nick site.
    pub rt_template: Rna,
    /// Complete 3'-extended pegRNA sequence (`spacer` + scaffold proxy + `rt_template` + `pbs`).
    pub pam_index: usize,
    /// Distance in nucleotides from the H840A nick site (`pam_index - 3`) to the edit position.
    pub nick_to_edit_distance: usize,
    /// Predicted Prime Editing efficiency ($0.0$ to $1.0$).
    pub estimated_efficiency: f64,
}

fn dna_to_reverse_complement_rna(dna_slice: &[u8]) -> Result<Rna> {
    let mut rna_bytes = Vec::with_capacity(dna_slice.len());
    for &b in dna_slice.iter().rev() {
        let comp = match b {
            b'A' => b'U',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        };
        rna_bytes.push(comp);
    }
    Rna::new(&rna_bytes)
}

/// Designs a Prime Editing (`PE2`) `pegRNA` to replace the base at `edit_index` in `locus_dna`
/// with `new_base`, using an SpCas9 `NGG` PAM located upstream of `edit_index` so that the H840A
/// nick site (`pam_index - 3`) lies 5' of `edit_index`.
pub fn design_prime_editing_pegrna(
    locus_dna: &Dna,
    edit_index: usize,
    new_base: u8,
    pbs_length: usize,
    homology_flap_length: usize,
) -> Result<Vec<PegRnaDesign>> {
    let bytes = &locus_dna.sequence;
    if edit_index >= bytes.len() {
        return Err(BioError::InvalidInput(format!(
            "edit_index ({edit_index}) is out of bounds for locus length ({})",
            bytes.len()
        )));
    }
    if !(8..=17).contains(&pbs_length) {
        return Err(BioError::InvalidInput(
            "Prime editing PBS length must be between 8 and 17 nt (optimal 10-13 nt)".to_string(),
        ));
    }

    let pams = find_spcas9_pam(locus_dna);
    let mut designs = Vec::new();

    for pam in pams {
        if pam.index < 20 {
            continue;
        }
        let nick_pos = pam.index - 3;
        // Edit must lie at or 3' of the nick site on the non-target strand, within 30 nt
        if edit_index < nick_pos || edit_index - nick_pos > 30 {
            continue;
        }
        if nick_pos < pbs_length || edit_index + 1 + homology_flap_length > bytes.len() {
            continue;
        }

        let proto_slice = &bytes[(pam.index - 20)..pam.index];
        let spacer_bytes: Vec<u8> = proto_slice
            .iter()
            .map(|&b| if b == b'T' { b'U' } else { b })
            .collect();
        let spacer = Rna::new(&spacer_bytes)?;

        // PBS hybridizes to the nicked strand [nick_pos - pbs_length .. nick_pos]
        let pbs_dna = &bytes[(nick_pos - pbs_length)..nick_pos];
        let pbs = dna_to_reverse_complement_rna(pbs_dna)?;

        // Edited DNA segment [nick_pos .. edit_index + 1 + homology_flap_length]
        let rtt_end = edit_index + 1 + homology_flap_length;
        let mut edited_dna_flap = bytes[nick_pos..rtt_end].to_vec();
        edited_dna_flap[edit_index - nick_pos] = new_base.to_ascii_uppercase();
        let rt_template = dna_to_reverse_complement_rna(&edited_dna_flap)?;

        let nick_to_edit = edit_index - nick_pos;
        let pbs_gc = pbs_dna
            .iter()
            .filter(|&&b| b == b'G' || b == b'C')
            .count() as f64
            / pbs_length as f64;
        let gc_factor = (1.0 - (pbs_gc - 0.55).abs() * 1.2).clamp(0.4, 1.0);
        let dist_factor = (1.0 - (nick_to_edit as f64 / 35.0) * 0.55).clamp(0.3, 1.0);
        let disrupts_pam_or_seed = (pam.index..=pam.index + 2).contains(&edit_index)
            || (pam.index - 6..pam.index).contains(&edit_index);
        let pam_disruption_boost = if disrupts_pam_or_seed { 1.15 } else { 0.90 };
        let estimated_efficiency =
            (0.65 * gc_factor * dist_factor * pam_disruption_boost).clamp(0.10, 0.88);

        designs.push(PegRnaDesign {
            spacer,
            pbs,
            rt_template,
            pam_index: pam.index,
            nick_to_edit_distance: nick_to_edit,
            estimated_efficiency,
        });
    }

    designs.sort_by(|a, b| {
        b.estimated_efficiency
            .partial_cmp(&a.estimated_efficiency)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(designs)
}

/// Evaluates chromosomal translocation probability ($0.0$ to $1.0$) during multiplex genome editing.
///
/// Simultaneous double-strand breaks ($N_{\text{DSB}}$) generate $\binom{N_{\text{DSB}}}{2}$
/// pairwise translocation opportunities, whereas DSB-free Base/Prime nickase edits ($N_{\text{nick}}$)
/// incur $>90\%$ lower translocation risk (Bothmer et al. 2017 *Nat Commun* 8:13905).
pub fn evaluate_multiplex_translocation_risk(
    simultaneous_dsb_edits: usize,
    simultaneous_nickase_edits: usize,
) -> f64 {
    let dsb_pairs = if simultaneous_dsb_edits >= 2 {
        (simultaneous_dsb_edits * (simultaneous_dsb_edits - 1)) / 2
    } else {
        0
    };
    let nick_pairs = if simultaneous_nickase_edits >= 2 {
        (simultaneous_nickase_edits * (simultaneous_nickase_edits - 1)) / 2
    } else {
        0
    };

    // Per-pair inter-chromosomal translocation probability: ~1.8% per simultaneous DSB pair vs ~0.08% per nickase pair
    let raw_prob = 1.0 - (0.982_f64.powi(dsb_pairs as i32) * 0.9992_f64.powi(nick_pairs as i32));
    raw_prob.clamp(0.0, 0.99)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_abe8e_base_editor_design_and_bystander_penalty() {
        // 23-nt sequence with PAM "AGG" at index 20:
        // Protospacer positions 1..=20 (0-based 0..20)
        // Place a single 'A' at 1-based position 6 (0-based index 5), surrounded by C/G/T
        let clean_locus = Dna::new(b"CCCGCTACCCTGCTCCCTGCAGG").unwrap();
        let clean_hits = design_base_editor_guides(&clean_locus, 6, BaseEditorType::Abe8e).unwrap();
        assert_eq!(clean_hits.len(), 1);
        assert_eq!(clean_hits[0].target_protospacer_pos_1idx, 7);
        assert_eq!(clean_hits[0].bystander_count, 0);
        assert!(clean_hits[0].precision_score > 0.80);

        // Now add 2 bystander 'A's at indices 4 and 5 (1-based positions 5 and 6)
        let bystander_locus = Dna::new(b"CCCGAAACCCTGCTCCCTGCAGG").unwrap();
        let bystander_hits =
            design_base_editor_guides(&bystander_locus, 6, BaseEditorType::Abe8e).unwrap();
        assert_eq!(bystander_hits.len(), 1);
        assert_eq!(bystander_hits[0].bystander_count, 2);
        assert!(
            bystander_hits[0].precision_score < clean_hits[0].precision_score * 0.5,
            "Bystander A nucleotides in window must penalize precision score"
        );
    }

    #[test]
    fn test_prime_editing_pegrna_design() {
        // PAM "GGG" at index 20 -> H840A nick at index 17
        // Edit target at index 19 (2 nt downstream of nick), plus 12-nt downstream sequence
        let locus = Dna::new(b"ATGCGATCGATCGATCGATCGGGATCGATCGATCGA").unwrap();
        let designs = design_prime_editing_pegrna(&locus, 19, b'T', 12, 10).unwrap();
        assert!(!designs.is_empty());
        let best = &designs[0];
        assert_eq!(best.spacer.sequence.len(), 20);
        assert_eq!(best.pbs.sequence.len(), 12);
        assert_eq!(best.nick_to_edit_distance, 2);
        assert_eq!(best.rt_template.sequence.len(), 13); // (19 - 17) + 1 + 10 = 13 nt
        assert!(best.estimated_efficiency > 0.40);
    }

    #[test]
    fn test_multiplex_dsb_vs_nickase_translocation_risk() {
        // 6 simultaneous Cas9 DSBs -> 15 pairwise translocation combinations (~24% translocation risk)
        let risk_6_dsb = evaluate_multiplex_translocation_risk(6, 0);
        // 6 simultaneous Base/Prime nickase edits -> >90% lower translocation risk
        let risk_6_nick = evaluate_multiplex_translocation_risk(0, 6);
        assert!(risk_6_dsb > 0.20);
        assert!(risk_6_nick < 0.02);
        assert!(risk_6_nick < risk_6_dsb * 0.10);
    }
}
