//! Example: End-to-end Longevity & Synthetic Immortal Genome Design from Real Public Genome Data.
//!
//! Demonstrates:
//! 1. Loading actual terminal telomere and promoter sequences from T2T-CHM13 assembly files.
//! 2. Profiling canonical TTAGGG / CCCTAA repeats and subtelomeric variant transitions.
//! 3. Running Design Engine 2.0 (High-Fidelity) with exact coordinates, PAMs, and CFD evaluations.
//! 4. Enforcing clinical oncogenic safety locks and sequential intervention ordering.

use scies_bio_th::aging::{
    design_high_fidelity_longevity_package, load_telomere_regions_from_file,
    HighFidelityAgeParams, LongevityStrategy,
};
use scies_bio_th::io::read_fasta_record_matching;
use scies_bio_th::seq::Dna;
use std::fs::File;
use std::io::BufReader;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🧬 scies-bio-th — Real Genome Integration & High-Fidelity Longevity Design (v0.5.0)\n");

    // 1. Load and analyze real terminal telomeres from T2T-CHM13 public genome sample
    let telomere_file = "testdata/t2t_chr5_telomere_tert.fasta";
    println!("1. Loading Chromosome 5 Terminal Region from '{}'...", telomere_file);
    let regions = load_telomere_regions_from_file(telomere_file, &["chr5"])?;

    if let Some(terminal) = regions.first() {
        println!("   Chromosome: {} | Arm: {:?}", terminal.chrom, terminal.arm);
        println!("   Scanned Length: {} bp | Strand: {:?}", terminal.sequence.sequence.len(), terminal.strand);
        println!("   Canonical Repeats: {} ({} bp)",
            terminal.analysis.canonical_repeats_count,
            terminal.analysis.canonical_repeat_len_bp
        );
        println!("   Subtelomeric Variant Repeats: {}", terminal.analysis.variant_repeats_count);
        println!("   Telomeric Density: {:.1}%", terminal.analysis.telomeric_density * 100.0);
        println!("   Subtelomeric Boundary Transition @ {} bp", terminal.analysis.subtelomeric_transition_bp);
        println!("   Shelterin Capping Index: {:.2}", terminal.analysis.estimated_shelterin_capping);
        println!("   Critically Short Hayflick Alarm: {}", terminal.analysis.is_critically_short);
    }
    println!();

    // 2. Extract Real Loci Sequences for Precision sgRNA Targeting
    println!("2. Extracting Real Genomic Target Loci from T2T-CHM13...");
    
    // TERT promoter from chr5 FASTA
    let f1 = File::open("testdata/t2t_chr5_telomere_tert.fasta")?;
    let r1 = BufReader::new(f1);
    let tert_record = read_fasta_record_matching(r1, |id| id.contains("tert_promoter"))?
        .ok_or("Failed to locate TERT promoter in testdata")?;
    let tert_dna = Dna::new(&tert_record.sequence)?;
    println!("   Extracted TERT Promoter locus: {} bp", tert_dna.sequence.len());

    // SIRT6 locus from chr19 FASTA
    let f2 = File::open("testdata/t2t_chr19_sirt6.fasta")?;
    let r2 = BufReader::new(f2);
    let sirt6_record = read_fasta_record_matching(r2, |id| id.contains("sirt6_locus"))?
        .ok_or("Failed to locate SIRT6 locus in testdata")?;
    let sirt6_dna = Dna::new(&sirt6_record.sequence)?;
    println!("   Extracted SIRT6 Regulatory locus: {} bp\n", sirt6_dna.sequence.len());

    // 3. Clinical Subject Baseline Parameters (Middle-aged subject with accelerated aging)
    let subject_baseline = HighFidelityAgeParams {
        chronological_age: 56.0,
        cpg_methylation_fraction: 0.54,
        telomere_length_bp: 5_900.0,
        gamma_h2ax_foci: 2.2, // Elevated DSB damage
        sasp_burden: 0.45,    // Chronic senescent secretome burden
    };

    println!("3. Running High-Fidelity Longevity Design Pipeline (Strategy: Balanced)...");
    let real_loci = [
        ("TERT", &tert_dna),
        ("SIRT6", &sirt6_dna),
    ];

    let report = design_high_fidelity_longevity_package(
        &subject_baseline,
        LongevityStrategy::Balanced,
        Some(&real_loci),
    );

    // 4. Output Detailed Structured Report
    println!("{}", report.export_detailed_summary());

    println!("✅ Real-genome high-fidelity longevity pipeline executed successfully!");
    Ok(())
}
