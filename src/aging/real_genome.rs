//! Real genome loading, telomere terminal sequence analysis, and subtelomeric boundary mapping.
//!
//! Provides direct support for extracting and profiling terminal chromosome regions from
//! public human genome assemblies (T2T-CHM13v2.0, GRCh38) without dummy approximations.

use crate::genome::Strand;
use crate::io::{FastaReader, FastaRecord};
use crate::seq::sha256::compute_sha256_hex;
use crate::seq::Dna;
use std::io::BufRead;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Canonical vertebrate telomere repeat units.
pub const CANONICAL_FORWARD_REPEAT: &[u8] = b"TTAGGG";
pub const CANONICAL_REVERSE_REPEAT: &[u8] = b"CCCTAA";

/// Frequent human subtelomeric variant repeats (T2T-CHM13 / GRCh38).
pub const FORWARD_VARIANT_REPEATS: &[&[u8]] = &[b"TGAGGG", b"TCAGGG", b"TTGGGG", b"TTAGGA"];
pub const REVERSE_VARIANT_REPEATS: &[&[u8]] = &[b"CCCTCA", b"CCCTGA", b"CCCCAA", b"TCCTAA"];

/// Chromosome arm designating the terminal locus orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ChromosomeArm {
    /// Short arm terminal (5' end of reference chromosome, predominantly `CCCTAA` forward).
    PArm,
    /// Long arm terminal (3' end of reference chromosome, predominantly `TTAGGG` forward).
    QArm,
    /// Undetermined or unanchored contig.
    Unknown,
}

/// Detailed profile of a chromosome terminal telomeric tract.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TelomereAnalysis {
    /// Count of canonical hexanucleotide repeats (`TTAGGG` or `CCCTAA`).
    pub canonical_repeats_count: usize,
    /// Total length in base pairs spanned by canonical repeats (`count * 6`).
    pub canonical_repeat_len_bp: usize,
    /// Count of degenerate subtelomeric variant repeats.
    pub variant_repeats_count: usize,
    /// 0-based offset where telomeric repeat density drops below threshold into subtelomere.
    pub subtelomeric_transition_bp: usize,
    /// Proportion of scanned sequence composed of telomeric or variant motifs.
    pub telomeric_density: f64,
    /// Estimated shelterin protein complex capping integrity index (0.0 to 1.0).
    pub estimated_shelterin_capping: f64,
    /// Whether the measured length falls below critical Hayflick limit threshold (< 3,800 bp).
    pub is_critically_short: bool,
}

/// A parsed terminal chromosome telomere region loaded from public genome FASTA.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct TelomereRegion {
    /// Chromosome identifier (e.g. "chr5", "chr19").
    pub chrom: String,
    /// Chromosome arm (p-arm or q-arm).
    pub arm: ChromosomeArm,
    /// Genomic start coordinate (0-based inclusive).
    pub start: u64,
    /// Genomic end coordinate (0-based exclusive).
    pub end: u64,
    /// Strand orientation.
    pub strand: Strand,
    /// Extracted DNA sequence.
    pub sequence: Dna,
    /// Comprehensive terminal repeat and capping analysis.
    pub analysis: TelomereAnalysis,
}

/// Analyzes real telomeric sequence data, detecting canonical & variant repeats and subtelomeric transitions.
pub fn analyze_real_telomere(sequence: &Dna) -> TelomereAnalysis {
    let bytes = &sequence.sequence;
    if bytes.len() < 6 {
        return TelomereAnalysis {
            canonical_repeats_count: 0,
            canonical_repeat_len_bp: 0,
            variant_repeats_count: 0,
            subtelomeric_transition_bp: 0,
            telomeric_density: 0.0,
            estimated_shelterin_capping: 0.1,
            is_critically_short: true,
        };
    }

    let mut forward_canonical = 0usize;
    let mut reverse_canonical = 0usize;
    let mut forward_variants = 0usize;
    let mut reverse_variants = 0usize;

    // 1. Phase-tolerant hexamer scan (handles arbitrary 5' coordinate offsets)
    let mut i = 0;
    while i + 6 <= bytes.len() {
        let hex = &bytes[i..i + 6];
        if hex == CANONICAL_FORWARD_REPEAT {
            forward_canonical += 1;
            i += 6;
        } else if hex == CANONICAL_REVERSE_REPEAT {
            reverse_canonical += 1;
            i += 6;
        } else if FORWARD_VARIANT_REPEATS.contains(&hex) {
            forward_variants += 1;
            i += 6;
        } else if REVERSE_VARIANT_REPEATS.contains(&hex) {
            reverse_variants += 1;
            i += 6;
        } else {
            i += 1;
        }
    }

    let (canonical_count, variant_count) = if reverse_canonical >= forward_canonical {
        (reverse_canonical, reverse_variants)
    } else {
        (forward_canonical, forward_variants)
    };

    let canonical_len = canonical_count * 6;
    let total_telomeric_bp = (canonical_count + variant_count) * 6;
    let telomeric_density = (total_telomeric_bp as f64 / bytes.len() as f64).min(1.0);

    // 2. Locate subtelomeric transition via phase-tolerant sliding windows.
    // Tracks the last window where telomeric repeat density >= 0.35 so that brief 5'-cap
    // degeneracy at the extreme chromosome tip (e.g. T2T-CHM13 chr5:1-200) does not mask
    // a 2,000+ bp canonical telomere tract immediately following it.
    let mut subtelomeric_transition = bytes.len();
    let window_size = 30.min(bytes.len());
    if window_size >= 12 {
        let step = 6;
        let mut pos = 0;
        let mut window_densities = Vec::new();
        while pos + window_size <= bytes.len() {
            let win = &bytes[pos..pos + window_size];
            let mut win_telomeric = 0;
            let mut j = 0;
            while j + 6 <= win.len() {
                let h = &win[j..j + 6];
                if h == CANONICAL_FORWARD_REPEAT
                    || h == CANONICAL_REVERSE_REPEAT
                    || FORWARD_VARIANT_REPEATS.contains(&h)
                    || REVERSE_VARIANT_REPEATS.contains(&h)
                {
                    win_telomeric += 6;
                    j += 6;
                } else {
                    j += 1;
                }
            }
            let density = win_telomeric as f64 / win.len() as f64;
            window_densities.push((pos, density));
            pos += step;
        }

        // Find the first low-density window after which telomeric density never recovers >= 0.50
        for idx in 0..window_densities.len() {
            if window_densities[idx].1 < 0.35 {
                let recovers_later = window_densities[idx + 1..]
                    .iter()
                    .any(|&(_, d)| d >= 0.50);
                if !recovers_later {
                    subtelomeric_transition = window_densities[idx].0;
                    break;
                }
            }
        }
    }

    let is_critically_short = canonical_len < 3_800;
    let shelterin_capping = if is_critically_short {
        ((canonical_len as f64 / 3_800.0) * 0.5).clamp(0.05, 0.50)
    } else {
        let len_score = ((canonical_len - 3_800) as f64 / 8_000.0).clamp(0.0, 1.0);
        let purity_score = (canonical_count as f64 / (canonical_count + variant_count).max(1) as f64).clamp(0.0, 1.0);
        (0.6 * len_score + 0.4 * purity_score).clamp(0.1, 0.99)
    };

    TelomereAnalysis {
        canonical_repeats_count: canonical_count,
        canonical_repeat_len_bp: canonical_len,
        variant_repeats_count: variant_count,
        subtelomeric_transition_bp: subtelomeric_transition,
        telomeric_density,
        estimated_shelterin_capping: shelterin_capping,
        is_critically_short,
    }
}

/// Parses FASTA records from a buffered reader and identifies all terminal telomere regions.
pub fn load_telomere_regions_from_fasta<R: BufRead>(reader: R) -> crate::Result<Vec<TelomereRegion>> {
    let fasta = FastaReader::new(reader)?;
    let mut regions = Vec::new();

    for res in fasta {
        let record: FastaRecord = res?;
        let id_lower = record.id.to_lowercase();
        
        let arm = if id_lower.contains("p_arm") || id_lower.contains("_p") || id_lower.contains("5prime") {
            ChromosomeArm::PArm
        } else if id_lower.contains("q_arm") || id_lower.contains("_q") || id_lower.contains("3prime") {
            ChromosomeArm::QArm
        } else {
            ChromosomeArm::Unknown
        };

        let strand = if arm == ChromosomeArm::PArm {
            Strand::Reverse
        } else {
            Strand::Forward
        };

        let seq_dna = Dna::new(&record.sequence)?;
        let analysis = analyze_real_telomere(&seq_dna);

        let end_coord = record.sequence.len() as u64;
        let chrom = record.id.split_whitespace().next().unwrap_or("unknown").to_string();

        regions.push(TelomereRegion {
            chrom,
            arm,
            start: 0,
            end: end_coord,
            strand,
            sequence: seq_dna,
            analysis,
        });
    }

    Ok(regions)
}

/// Convenience helper to load telomere regions from a FASTA file on disk, optionally filtering by chromosomes.
pub fn load_telomere_regions_from_file(
    file_path: &str,
    target_chromosomes: &[&str],
) -> crate::Result<Vec<TelomereRegion>> {
    let file = std::fs::File::open(file_path)
        .map_err(|e| crate::BioError::IoError(format!("Failed to open FASTA file '{}': {}", file_path, e)))?;
    let reader = std::io::BufReader::new(file);
    let all_regions = load_telomere_regions_from_fasta(reader)?;

    if target_chromosomes.is_empty() {
        return Ok(all_regions);
    }

    let filtered = all_regions
        .into_iter()
        .filter(|r| {
            target_chromosomes
                .iter()
                .any(|c| r.chrom.eq_ignore_ascii_case(c) || r.chrom.contains(c))
        })
        .collect();

    Ok(filtered)
}

/// Loads all genomic locus records from a buffered FASTA stream as `(identifier, Dna)`.
pub fn load_loci_from_fasta<R: BufRead>(reader: R) -> crate::Result<Vec<(String, Dna)>> {
    let fasta = FastaReader::new(reader)?;
    let mut loci = Vec::new();
    for res in fasta {
        let record = res?;
        let dna = Dna::new(&record.sequence)?;
        loci.push((record.id, dna));
    }
    Ok(loci)
}

/// Loads all genomic locus records from a FASTA file on disk.
pub fn load_loci_from_file(file_path: &str) -> crate::Result<Vec<(String, Dna)>> {
    let file = std::fs::File::open(file_path)
        .map_err(|e| crate::BioError::IoError(format!("Failed to open FASTA file '{}': {}", file_path, e)))?;
    let reader = std::io::BufReader::new(file);
    load_loci_from_fasta(reader)
}

/// Provenance and coordinate metadata associated with a real genomic locus extracted from reference FASTA.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct RealLocusMetadata {
    /// Genome assembly version (e.g., "T2T-CHM13v2.0" or "GRCh38.p14").
    pub assembly: String,
    /// Chromosome identifier (e.g., "chr5", "chr19").
    pub chrom: String,
    /// 0-based inclusive start coordinate.
    pub start: u64,
    /// 0-based exclusive end coordinate.
    pub end: u64,
    /// Strand orientation.
    pub strand: Strand,
    /// Cryptographic SHA-256 hash of the source file from which this locus was loaded.
    pub source_file_hash: String,
    /// HGNC gene symbol if annotated (e.g., Some("TERT")).
    pub gene_symbol: Option<String>,
}

/// A parsed real genomic locus with validated sequence and reference metadata.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct RealGenomicLocus {
    pub id: String,
    pub sequence: Dna,
    pub metadata: RealLocusMetadata,
}

/// Loads all real genomic locus records from a FASTA file on disk, extracting coordinate/assembly metadata
/// and recording the cryptographic SHA-256 hash of the source file.
pub fn load_real_genomic_loci_from_file(file_path: &str) -> crate::Result<Vec<RealGenomicLocus>> {
    let raw_bytes = std::fs::read(file_path)
        .map_err(|e| crate::BioError::IoError(format!("Failed to read FASTA file '{}': {}", file_path, e)))?;
    let file_hash = compute_sha256_hex(&raw_bytes);

    let reader = std::io::BufReader::new(raw_bytes.as_slice());
    let fasta = FastaReader::new(reader)?;
    let mut loci = Vec::new();

    for res in fasta {
        let record = res?;
        let dna = Dna::new(&record.sequence)?;
        let header = &record.id;

        // Parse assembly, e.g. [assembly=GRCh38.p14] or [assembly=T2T-CHM13v2.0]
        let assembly = if let Some(idx) = header.find("assembly=") {
            let rem = &header[idx + 9..];
            let end_idx = rem.find(']').unwrap_or(rem.len());
            rem[..end_idx].trim().to_string()
        } else if header.contains("T2T") {
            "T2T-CHM13v2.0".to_string()
        } else {
            "GRCh38.p14".to_string()
        };

        // Parse coordinates: chrX:start-end
        let (chrom, start, end) = parse_chrom_coordinates(&header).unwrap_or_else(|| {
            let chr = if record.id.contains("chr") {
                let parts: Vec<&str> = record.id.split('_').collect();
                parts.iter().find(|p| p.starts_with("chr")).copied().unwrap_or("chrUnknown").to_string()
            } else {
                "chrUnknown".to_string()
            };
            (chr, 1, 1 + dna.sequence.len() as u64)
        });

        // Parse gene symbol: [gene=...] or infer from id
        let gene_symbol = if let Some(idx) = header.find("gene=") {
            let rem = &header[idx + 5..];
            let end_idx = rem.find(']').unwrap_or(rem.len());
            let sym = rem[..end_idx].trim();
            Some(sym.to_string())
        } else {
            infer_gene_from_id(&record.id)
        };

        let metadata = RealLocusMetadata {
            assembly,
            chrom,
            start,
            end,
            strand: Strand::Forward,
            source_file_hash: file_hash.clone(),
            gene_symbol,
        };

        loci.push(RealGenomicLocus {
            id: record.id,
            sequence: dna,
            metadata,
        });
    }

    Ok(loci)
}

fn parse_chrom_coordinates(header: &str) -> Option<(String, u64, u64)> {
    for token in header.split_whitespace() {
        if let Some(colon_idx) = token.find(':') {
            let chrom = token[..colon_idx].trim_start_matches('>').to_string();
            let coords = &token[colon_idx + 1..];
            if let Some(dash_idx) = coords.find('-') {
                let start_str = coords[..dash_idx].replace(',', "");
                let end_str = coords[dash_idx + 1..].replace(',', "");
                if let (Ok(start), Ok(end)) = (start_str.parse::<u64>(), end_str.parse::<u64>()) {
                    return Some((chrom, start, end));
                }
            }
        }
    }
    None
}

fn infer_gene_from_id(id: &str) -> Option<String> {
    let lower = id.to_ascii_lowercase();
    for gene in &[
        "tert", "sirt6", "tp53", "cdkn2a", "foxo3", "parp1", "pot1", "lmna", "atm", "terc", "mtor",
    ] {
        if lower.contains(gene) {
            return Some(gene.to_ascii_uppercase());
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_analyze_real_telomere_canonical_and_variants() {
        let mut seq_bytes = Vec::new();
        // 10 canonical CCCTAA repeats (60 bp)
        for _ in 0..10 {
            seq_bytes.extend_from_slice(b"CCCTAA");
        }
        // 2 variant repeats (12 bp)
        seq_bytes.extend_from_slice(b"CCCTGA");
        seq_bytes.extend_from_slice(b"CCCCAA");
        // Non-telomeric flanking sequence (30 bp)
        seq_bytes.extend_from_slice(b"GATCGATCGATCGATCGATCGATCGATCGA");

        let dna = Dna::new(&seq_bytes).unwrap();
        let analysis = analyze_real_telomere(&dna);

        assert_eq!(analysis.canonical_repeats_count, 10);
        assert_eq!(analysis.canonical_repeat_len_bp, 60);
        assert_eq!(analysis.variant_repeats_count, 2);
        assert!(analysis.telomeric_density > 0.65);
        assert!(analysis.subtelomeric_transition_bp <= 72);
    }

    #[test]
    fn test_load_telomere_regions_from_fasta_stream() {
        let fasta_content = ">chr5_telomere_p_arm chr5:1-120\n\
                             CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA\n\
                             CCCTGACCCCAAAGTCGATCGATCGATCGATCGATCGATCGATCGATC\n";
        let reader = std::io::BufReader::new(fasta_content.as_bytes());
        let regions = load_telomere_regions_from_fasta(reader).unwrap();

        assert_eq!(regions.len(), 1);
        let reg = &regions[0];
        assert_eq!(reg.arm, ChromosomeArm::PArm);
        assert_eq!(reg.strand, Strand::Reverse);
        assert_eq!(reg.analysis.canonical_repeats_count, 10);
        assert_eq!(reg.analysis.variant_repeats_count, 2);
    }

    #[test]
    fn test_load_real_genomic_loci_with_metadata_and_provenance_hash() {
        let loci = load_real_genomic_loci_from_file("testdata/t2t_chr5_telomere_tert.fasta").unwrap();
        assert!(!loci.is_empty());
        let tert = loci.iter().find(|l| l.id.contains("tert")).unwrap();
        assert_eq!(tert.metadata.chrom, "chr5");
        assert_eq!(tert.metadata.assembly, "GRCh38.p14");
        assert_eq!(tert.metadata.start, 1294801);
        assert_eq!(tert.metadata.end, 1295600);
        assert_eq!(tert.metadata.gene_symbol.as_deref(), Some("ENSG00000164362"));
        assert_eq!(tert.metadata.source_file_hash.len(), 64);
    }
}
