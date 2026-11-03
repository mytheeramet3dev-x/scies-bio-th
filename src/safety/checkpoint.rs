//! p53 / Rb / ATM DNA damage and replication-stress checkpoint activation model.
//!
//! Models ATM/ATR kinase phosphorylation of p53 (Ser15/Ser20), p21-mediated CDK2/Rb G1/S
//! cell-cycle arrest, and PUMA/BAX apoptotic threshold commitment.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Upstream DNA damage, telomere uncapping, and oncogenic stress signals feeding into ATM/p53/Rb.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CheckpointInput {
    /// Double-strand break burden ($\gamma$-H2AX foci per cell, typically $0.1\text{--}15.0$).
    pub dsb_foci_per_cell: f64,
    /// Shortest telomere length in base pairs (critical shelterin uncapping below $4,000\text{ bp}$).
    pub shortest_telomere_bp: f64,
    /// Relative p53 (TP53) expression multiplier ($1.0 = \text{wild-type diploid}$).
    pub tp53_expression_fold: f64,
    /// Relative p16INK4a / Rb pathway integrity ($1.0 = \text{intact}$).
    pub rb_p16_integrity: f64,
    /// Mitogenic / proliferative drive fold (e.g. MYC/RAS or high TERT + OSK).
    pub mitogenic_drive_fold: f64,
}

/// Cellular fate dictated by p53 / Rb checkpoint integration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum CellFateDecision {
    /// Normal cell-cycle progression (low damage, intact checkpoints).
    ProliferationPermitted,
    /// Transient G1/S arrest via p53 -> p21 -> Rb hypophosphorylation for DNA repair.
    TransientG1ArrestRepair,
    /// Permanent replicative or stress-induced senescence.
    SenescenceLock,
    /// p53/PUMA/BAX-mediated apoptosis (safe elimination of damaged clone).
    ApoptoticClearance,
    /// Checkpoint bypass crisis (high damage or mitogenic drive with deficient p53/Rb).
    GenomicInstabilityBypass,
}

/// Quantitative activation state of the ATM/p53/Rb checkpoint network.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CheckpointResponse {
    /// ATM/ATR sensor kinase activation score in $[0.0, 1.0]$.
    pub atm_atr_activity: f64,
    /// Stabilized p53 transcriptional activity score in $[0.0, 1.0]$.
    pub p53_activation_score: f64,
    /// Active (hypophosphorylated) Rb G1-arrest barrier in $[0.0, 1.0]$.
    pub rb_arrest_barrier: f64,
    /// Threshold above which p53 triggers cell-cycle arrest ($0.35$).
    pub arrest_threshold: f64,
    /// Resulting cellular fate decision.
    pub fate: CellFateDecision,
}

/// Evaluates ATM/p53/Rb pathway activation and cell-fate commitment.
pub fn evaluate_checkpoint_activation(input: &CheckpointInput) -> CheckpointResponse {
    let dsb_signal = (input.dsb_foci_per_cell.max(0.0) / 5.0).clamp(0.0, 1.0);
    let telomere_uncapping = if input.shortest_telomere_bp < 4500.0 {
        ((4500.0 - input.shortest_telomere_bp) / 2500.0).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let arf_oncogenic_signal = ((input.mitogenic_drive_fold - 1.2).max(0.0) / 2.5).clamp(0.0, 1.0);

    // ATM/ATR responds to DSBs and uncapped telomeres
    let atm_atr_activity = (0.65 * dsb_signal + 0.55 * telomere_uncapping).clamp(0.0, 1.0);

    // p53 stabilization integrates ATM/ATR phosphorylation + p14ARF (oncogenic stress) scaled by TP53 dosage
    let tp53_capacity = (input.tp53_expression_fold.max(0.0) / 1.5).clamp(0.0, 1.5);
    let raw_p53 = (0.75 * atm_atr_activity + 0.45 * arf_oncogenic_signal) * tp53_capacity;
    let p53_activation_score = raw_p53.clamp(0.0, 1.0);

    // Rb barrier integrates p21 (from p53) and p16INK4a against mitogenic cyclin D/E drive
    let p21_induction = 0.60 * p53_activation_score;
    let rb_arrest_barrier = ((p21_induction + 0.50 * input.rb_p16_integrity.clamp(0.0, 1.5))
        / input.mitogenic_drive_fold.max(0.5))
    .clamp(0.0, 1.0);

    let arrest_threshold = 0.35;
    let severe_stress = atm_atr_activity > 0.45 || arf_oncogenic_signal > 0.50;

    let fate = if severe_stress && input.tp53_expression_fold < 0.4 && rb_arrest_barrier < 0.30 {
        CellFateDecision::GenomicInstabilityBypass
    } else if p53_activation_score >= 0.75 {
        CellFateDecision::ApoptoticClearance
    } else if p53_activation_score >= arrest_threshold && telomere_uncapping > 0.5 {
        CellFateDecision::SenescenceLock
    } else if p53_activation_score >= arrest_threshold {
        CellFateDecision::TransientG1ArrestRepair
    } else {
        CellFateDecision::ProliferationPermitted
    };

    CheckpointResponse {
        atm_atr_activity,
        p53_activation_score,
        rb_arrest_barrier,
        arrest_threshold,
        fate,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_p53_checkpoint_activation_above_threshold_on_dna_damage() {
        // Validation target: p53 activation score > threshold on simulated DNA damage signal
        let damaged = CheckpointInput {
            dsb_foci_per_cell: 4.5,
            shortest_telomere_bp: 3800.0,
            tp53_expression_fold: 1.5,
            rb_p16_integrity: 1.0,
            mitogenic_drive_fold: 1.0,
        };
        let resp = evaluate_checkpoint_activation(&damaged);
        assert!(resp.p53_activation_score > resp.arrest_threshold);
        assert!(
            resp.fate == CellFateDecision::ApoptoticClearance
                || resp.fate == CellFateDecision::TransientG1ArrestRepair
        );
    }
}
