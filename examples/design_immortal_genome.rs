//! Example: Designing a Synthetic Immortal Genome Package with scies-bio-th.
//!
//! Run with: `cargo run --example design_immortal_genome`

use scies_bio_th::aging::{
    apply_telomerase, count_telomeric_repeats, design_immortal_genome_package,
    estimate_biological_age, predict_replicative_lifespan, simulate_telomere_shortening,
    AgingProfile, LongevityStrategy, TelomeraseActivity,
};
use scies_bio_th::seq::Dna;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🧬 scies-bio-th — Longevity & Synthetic Genome Design (v0.4.0)\n");

    // 1. Terminal Chromosome Telomere Scanning
    let sample_telomere_dna = Dna::new(b"ATGCGTTAGGGTTAGGGTTAGGGCCCTAAGATCGATC")?;
    let repeat_count = count_telomeric_repeats(&sample_telomere_dna);
    println!("1. Terminal Repeat Analysis:");
    println!("   Sequence Length: {} bp", sample_telomere_dna.sequence.len());
    println!("   Detected Telomeric Repeats (TTAGGG / CCCTAA): {}", repeat_count);

    // 2. Baseline Aging Profile Simulation
    let chronological_age = 52.0;
    // Simulate somatic telomere erosion after 65 cellular divisions with mild oxidative stress
    let baseline_telomere = simulate_telomere_shortening(12_000.0, 65, 0.4);
    let baseline_bio_age = estimate_biological_age(
        chronological_age,
        0.55, // Epigenetic methylation drift
        baseline_telomere.current_length_bp,
    );

    let baseline_profile = AgingProfile {
        chronological_age_years: chronological_age,
        biological_age_years: baseline_bio_age,
        telomere: baseline_telomere.clone(),
        senescence_burden: 0.35,
        dna_damage_foci: 1.8,
        risk_scores: std::collections::HashMap::new(),
    };

    println!("\n2. Baseline Subject Assessment:");
    println!("   Chronological Age: {:.1} years", baseline_profile.chronological_age_years);
    println!("   Biological Epigenetic Age: {:.1} years", baseline_profile.biological_age_years);
    println!("   Telomere Length: {:.0} bp", baseline_profile.telomere.current_length_bp);
    println!(
        "   Remaining Replicative Lifespan: {} divisions",
        predict_replicative_lifespan(&baseline_profile.telomere)
    );

    // 3. Synthetic Genome Design Pipeline (Balanced Strategy)
    println!("\n3. Running Immortal Genome Engineering Pipeline (Balanced Strategy)...");
    let report = design_immortal_genome_package(&baseline_profile, LongevityStrategy::Balanced);

    println!("\n{}", report.export_summary());

    // 4. In Silico Telomerase Rejuvenation Simulation
    println!("4. In Silico Telomere Maintenance Simulation:");
    let mut simulated_telomere = baseline_profile.telomere.clone();
    let telomerase = TelomeraseActivity::synthetic_enhanced(true);
    let added_bp = apply_telomerase(&mut simulated_telomere, &telomerase, 30);

    println!("   Telomerase Cycles Applied: 30");
    println!("   Net Nucleotides Added: +{:.0} bp", added_bp);
    println!("   Restored Telomere Length: {:.0} bp", simulated_telomere.current_length_bp);
    println!(
        "   Expanded Replicative Lifespan: {} divisions",
        predict_replicative_lifespan(&simulated_telomere)
    );
    println!("   End Protection Integrity: {:.1}%", simulated_telomere.end_protection_score() * 100.0);

    println!("\n✅ Longevity design pipeline completed successfully!");
    Ok(())
}
