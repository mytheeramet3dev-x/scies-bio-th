use crate::genome::Interval;

/// The biological type of a genomic feature.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FeatureType {
    Gene,
    Transcript,
    Exon,
    Cds,
    FivePrimeUtr,
    ThreePrimeUtr,
    Other(String),
}

impl FeatureType {
    pub fn from_str(s: &str) -> Self {
        match s {
            "gene" => FeatureType::Gene,
            "transcript" | "mRNA" => FeatureType::Transcript,
            "exon" => FeatureType::Exon,
            "CDS" => FeatureType::Cds,
            "five_prime_UTR" | "5UTR" => FeatureType::FivePrimeUtr,
            "three_prime_UTR" | "3UTR" => FeatureType::ThreePrimeUtr,
            other => FeatureType::Other(other.to_string()),
        }
    }
}

/// A single annotated feature on a genome.
#[derive(Debug, Clone)]
pub struct Feature {
    pub feature_type: FeatureType,
    pub interval: Interval,
    pub source: Option<String>,
    pub score: Option<f64>,
    pub id: Option<String>,
    pub parent: Option<String>,
}

impl Feature {
    pub fn new(
        feature_type: FeatureType,
        interval: Interval,
    ) -> Self {
        Self {
            feature_type,
            interval,
            source: None,
            score: None,
            id: None,
            parent: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::Strand;

    #[test]
    fn test_feature_creation() {
        let iv = Interval::new("chr1", 1000, 5000, Strand::Forward).unwrap();
        let feat = Feature::new(FeatureType::Gene, iv);
        assert_eq!(feat.feature_type, FeatureType::Gene);
        assert_eq!(feat.interval.len(), 4000);
    }

    #[test]
    fn test_feature_type_from_str() {
        assert_eq!(FeatureType::from_str("gene"), FeatureType::Gene);
        assert_eq!(FeatureType::from_str("exon"), FeatureType::Exon);
        assert_eq!(FeatureType::from_str("CDS"), FeatureType::Cds);
    }
}
