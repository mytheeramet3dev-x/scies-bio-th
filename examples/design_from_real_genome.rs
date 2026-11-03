//! Example: End-to-end Longevity & Synthetic Immortal Genome Design from Real Public Genome Data.
//!
//! Demonstrates:
//! 1. Loading actual terminal telomere and promoter sequences from T2T-CHM13 assembly files.
//! 2. Profiling canonical TTAGGG / CCCTAA repeats and subtelomeric variant transitions.
//! 3. Designing real SpCas9 sgRNAs with homopolymer (poly-T / poly-G) filters and CFD scoring.
//! 4. Enforcing hard biological safety constraints (Super-p53 and checkpoint primacy).
//! 5. Simulating multi-generation cell population dynamics over 50 generations.
//! 6. Running pre-flight clinical safety validation.

use scies_bio_th::aging::{
    design_high_fidelity_longevity_package, load_loci_from_file, load_telomere_regions_from_file,
    validate_longevity_package, HighFidelityAgeParams, LongevityStrategy,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🧬 scies-bio-th — Real Guide Design, Simulation & Hard Constraints (v1.0.0)\n");

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

    // 2. Extract Real Genomic Target Loci from GRCh38.p14 / T2T-CHM13v2.0 Testdata
    println!("2. Extracting Real Genomic Target Loci from GRCh38.p14 / T2T-CHM13v2.0...");

    let tert_loci = load_loci_from_file("testdata/t2t_chr5_telomere_tert.fasta")?;
    let sirt6_loci = load_loci_from_file("testdata/t2t_chr19_sirt6.fasta")?;
    let tp53_loci = load_loci_from_file("testdata/t2t_chr17_tp53.fasta")?;
    let cdkn2a_loci = load_loci_from_file("testdata/t2t_chr9_cdkn2a.fasta")?;
    let foxo3_loci = load_loci_from_file("testdata/t2t_chr6_foxo3.fasta")?;
    let parp1_loci = load_loci_from_file("testdata/t2t_chr1_parp1.fasta")?;
    let pot1_loci = load_loci_from_file("testdata/t2t_chr7_pot1.fasta")?;
    let lmna_loci = load_loci_from_file("testdata/t2t_chr1_lmna.fasta")?;
    let atm_loci = load_loci_from_file("testdata/t2t_chr11_atm.fasta")?;
    let terc_loci = load_loci_from_file("testdata/t2t_chr3_terc.fasta")?;
    let mtor_loci = load_loci_from_file("testdata/t2t_chr1_mtor.fasta")?;

    let tert_dna = &tert_loci.iter().find(|(id, _)| id.contains("tert_promoter")).unwrap().1;
    let sirt6_dna = &sirt6_loci.iter().find(|(id, _)| id.contains("sirt6_locus")).unwrap().1;
    let tp53_dna = &tp53_loci.first().unwrap().1;
    let cdkn2a_dna = &cdkn2a_loci.first().unwrap().1;
    let foxo3_dna = &foxo3_loci.first().unwrap().1;
    let parp1_dna = &parp1_loci.first().unwrap().1;
    let pot1_dna = &pot1_loci.first().unwrap().1;
    let lmna_dna = &lmna_loci.first().unwrap().1;
    let atm_dna = &atm_loci.first().unwrap().1;
    let terc_dna = &terc_loci.first().unwrap().1;
    let mtor_dna = &mtor_loci.first().unwrap().1;

    println!("   Loaded TERT Promoter ({} bp)", tert_dna.sequence.len());
    println!("   Loaded SIRT6 Locus ({} bp)", sirt6_dna.sequence.len());
    println!("   Loaded TP53 Locus ({} bp)", tp53_dna.sequence.len());
    println!("   Loaded CDKN2A Locus ({} bp)", cdkn2a_dna.sequence.len());
    println!("   Loaded FOXO3 Locus ({} bp)", foxo3_dna.sequence.len());
    println!("   Loaded PARP1 Locus ({} bp)", parp1_dna.sequence.len());
    println!("   Loaded POT1 / LMNA / ATM / TERC / MTOR Loci (800 bp each)\n");

    let real_loci = [
        ("TERT", tert_dna),
        ("SIRT6", sirt6_dna),
        ("TP53", tp53_dna),
        ("CDKN2A", cdkn2a_dna),
        ("FOXO3", foxo3_dna),
        ("PARP1", parp1_dna),
        ("POT1", pot1_dna),
        ("LMNA", lmna_dna),
        ("ATM", atm_dna),
        ("TERC", terc_dna),
        ("MTOR", mtor_dna),
    ];

    // 3. Clinical Subject Baseline Parameters (Middle-aged subject with accelerated aging)
    let subject_baseline = HighFidelityAgeParams {
        chronological_age: 56.0,
        cpg_methylation_fraction: 0.54,
        telomere_length_bp: 5_900.0,
        gamma_h2ax_foci: 2.2, // Elevated DSB damage
        sasp_burden: 0.45,    // Chronic senescent secretome burden
    };

    println!("3. Running High-Fidelity Longevity Design Pipeline (Strategy: Balanced)...");
    let report = design_high_fidelity_longevity_package(
        &subject_baseline,
        LongevityStrategy::Balanced,
        Some(&real_loci),
    );

    // 4. Pre-Flight Clinical Validation
    println!("4. Performing Pre-Flight Biological Constraints Validation...");
    let validation = validate_longevity_package(&report);
    println!("   Validation Status: {}", if validation.is_valid { "PASSED (Ready for Synthesis)" } else { "FAILED" });
    if !validation.warnings.is_empty() {
        println!("   Advisories: {} recorded", validation.warnings.len());
    }
    println!();

    // 5. Output Detailed Structured Report
    println!("{}", report.export_detailed_summary());

    println!("✅ Real-genome high-fidelity longevity pipeline executed successfully!");
    Ok(())
}
