//! Telomere dynamics, end-replication problem simulation, and Shelterin complex capping models.
//!
//! Human chromosomes terminate with tandem arrays of the hexanucleotide repeat `5'-TTAGGG-3'`.
//! Somatic replication progressively erodes telomeres by ~50-100 bp per division, exacerbated
//! by oxidative stress, until reaching a critical Hayflick threshold that triggers replicative senescence.

use crate::seq::Dna;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Canonical vertebrate telomere repeat sequence (`5'-TTAGGG-3'`).
pub const VERTEBRATE_TELOMERE_REPEAT: &[u8] = b"TTAGGG";

/// Reverse complement telomeric repeat (`5'-CCCTAA-3'`).
pub const VERTEBRATE_TELOMERE_REVERSE: &[u8] = b"CCCTAA";

/// Current physiological and structural state of chromosomal telomeres.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TelomereState {
    /// Average terminal telomere length in base pairs (typical human neonate: ~12,000-15,000 bp).
    pub current_length_bp: f64,
    /// Basal shortening rate in bp per cellular division (typically ~50-70 bp).
    pub shortening_rate_bp: f64,
    /// Critical length threshold in bp triggering senescence / DDR alarm (typically ~3,500-4,000 bp).
    pub critical_length_bp: f64,
    /// Shelterin protein complex capping integrity index (0.0 to 1.0).
    pub shelterin_integrity: f64,
    /// Total cumulative cell divisions elapsed.
    pub divisions_elapsed: u32,
    /// Whether telomeres have breached the critical threshold into irreversible senescence.
    pub is_critically_short: bool,
}

impl TelomereState {
    /// Creates a baseline youthful telomere state (e.g. 12,000 bp).
    pub fn youthful_baseline() -> Self {
        Self {
            current_length_bp: 12_000.0,
            shortening_rate_bp: 60.0,
            critical_length_bp: 3_800.0,
            shelterin_integrity: 0.98,
            divisions_elapsed: 0,
            is_critically_short: false,
        }
    }

    /// Evaluates the end-capping protection against DNA damage response (DDR / ATM-ATR kinase activation).
    pub fn end_protection_score(&self) -> f64 {
        if self.is_critically_short {
            return 0.05;
        }
        let length_factor = ((self.current_length_bp - self.critical_length_bp) / 8_000.0).clamp(0.0, 1.0);
        (length_factor * 0.6 + self.shelterin_integrity * 0.4).clamp(0.0, 1.0)
    }
}

/// Simulates telomere attrition over a number of cellular divisions under oxidative stress.
///
/// `oxidative_stress`: 0.0 (minimal ROS) to 2.0+ (severe oxidative damage accelerating loss).
pub fn simulate_telomere_shortening(
    initial_length_bp: f64,
    divisions: u32,
    oxidative_stress: f64,
) -> TelomereState {
    let basal_rate = 60.0;
    let critical_length = 3_800.0;
    let stress_multiplier = 1.0 + oxidative_stress.clamp(0.0, 5.0) * 0.75;
    let effective_rate = basal_rate * stress_multiplier;

    let total_attrition = effective_rate * divisions as f64;
    let current_length = (initial_length_bp - total_attrition).max(1_000.0);
    let is_critically_short = current_length <= critical_length;

    // Shelterin uncapping accelerates once telomeres fall below 5,000 bp
    let shelterin_integrity = if current_length > 6_000.0 {
        0.98
    } else {
        ((current_length - 1_000.0) / 5_000.0).clamp(0.1, 0.95)
    };

    TelomereState {
        current_length_bp: current_length,
        shortening_rate_bp: effective_rate,
        critical_length_bp: critical_length,
        shelterin_integrity,
        divisions_elapsed: divisions,
        is_critically_short,
    }
}

/// Predicts remaining cellular replicative capacity (Hayflick limit) until critical arrest.
pub fn predict_replicative_lifespan(state: &TelomereState) -> u32 {
    if state.current_length_bp <= state.critical_length_bp || state.shortening_rate_bp <= 0.0 {
        return 0;
    }
    let remaining_bp = state.current_length_bp - state.critical_length_bp;
    (remaining_bp / state.shortening_rate_bp).floor() as u32
}

/// Scans a DNA sequence and counts exact occurrences of vertebrate telomeric repeats (`TTAGGG` and `CCCTAA`).
pub fn count_telomeric_repeats(dna: &Dna) -> usize {
    let bytes = &dna.sequence;
    if bytes.len() < 6 {
        return 0;
    }
    let mut count = 0;
    for i in 0..=(bytes.len() - 6) {
        let hexamer = &bytes[i..i + 6];
        if hexamer == VERTEBRATE_TELOMERE_REPEAT || hexamer == VERTEBRATE_TELOMERE_REVERSE {
            count += 1;
        }
    }
    count
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_telomere_shortening_and_hayflick_limit() {
        let state = simulate_telomere_shortening(10_000.0, 50, 0.0);
        // 10,000 - 50 * 60 = 7,000 bp
        assert_eq!(state.current_length_bp, 7_000.0);
        assert!(!state.is_critically_short);

        let remaining = predict_replicative_lifespan(&state);
        // (7,000 - 3,800) / 60 = 3,200 / 60 = 53 divisions
        assert_eq!(remaining, 53);
    }

    #[test]
    fn test_oxidative_stress_acceleration() {
        let calm = simulate_telomere_shortening(10_000.0, 40, 0.0);
        let stressed = simulate_telomere_shortening(10_000.0, 40, 1.0);
        assert!(stressed.current_length_bp < calm.current_length_bp);
    }

    #[test]
    fn test_count_telomeric_repeats() {
        let dna = Dna::new(b"ATGCGATTAGGGTTAGGGCCCTAAGATC").unwrap();
        // Has two TTAGGG and one CCCTAA = 3 telomeric repeats
        let c = count_telomeric_repeats(&dna);
        assert_eq!(c, 3);
    }
}
