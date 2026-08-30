use crate::genome::{Interval, Strand};
use crate::annotation::{Feature, FeatureType};
use crate::BioError;
use std::io::BufRead;

/// A single record from a GFF3 file.
#[derive(Debug, Clone)]
pub struct GffRecord {
    pub seqid: String,
    pub source: String,
    pub feature_type: String,
    pub start: u64,    // 1-based in GFF3, stored as 0-based internally
    pub end: u64,
    pub score: Option<f64>,
    pub strand: Strand,
    pub phase: Option<u8>,
    pub attributes: String,
}

impl GffRecord {
    /// Convert GFF3 record to annotation Feature.
    pub fn to_feature(&self) -> crate::Result<Feature> {
        let iv = Interval::new(&self.seqid, self.start, self.end, self.strand)?;
        let ftype = FeatureType::from_str(&self.feature_type);

        let mut feat = Feature::new(ftype, iv);
        feat.source = Some(self.source.clone());
        feat.score = self.score;

        // Simple attribute parsing for ID= and Parent=
        for attr in self.attributes.split(';') {
            let attr = attr.trim();
            if let Some(id) = attr.strip_prefix("ID=") {
                feat.id = Some(id.to_string());
            } else if let Some(parent) = attr.strip_prefix("Parent=") {
                feat.parent = Some(parent.to_string());
            }
        }

        Ok(feat)
    }
}

/// A streaming GFF3 parser.
pub struct GffReader<R> {
    reader: R,
}

impl<R: BufRead> GffReader<R> {
    pub fn new(reader: R) -> Self {
        Self { reader }
    }
}

fn parse_strand(s: &str) -> Strand {
    match s {
        "+" => Strand::Forward,
        "-" => Strand::Reverse,
        _ => Strand::Unknown,
    }
}

impl<R: BufRead> Iterator for GffReader<R> {
    type Item = crate::Result<GffRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line = String::new();
        loop {
            line.clear();
            match self.reader.read_line(&mut line) {
                Ok(0) => return None,
                Err(e) => return Some(Err(BioError::ParseError(e.to_string()))),
                _ => {}
            }
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue; // skip comments and blank lines
            }

            let fields: Vec<&str> = trimmed.splitn(9, '\t').collect();
            if fields.len() < 8 {
                return Some(Err(BioError::ParseError(
                    format!("GFF3: expected 9 fields, got {}", fields.len())
                )));
            }

            let start = match fields[3].parse::<u64>() {
                Ok(v) => v.saturating_sub(1), // GFF3 is 1-based, convert to 0-based
                Err(_) => return Some(Err(BioError::ParseError("Invalid start".into()))),
            };
            let end = match fields[4].parse::<u64>() {
                Ok(v) => v,
                Err(_) => return Some(Err(BioError::ParseError("Invalid end".into()))),
            };
            let score = match fields[5] {
                "." => None,
                s => s.parse::<f64>().ok(),
            };
            let phase = match fields[7] {
                "." => None,
                s => s.parse::<u8>().ok(),
            };
            let attributes = if fields.len() == 9 { fields[8].to_string() } else { String::new() };

            return Some(Ok(GffRecord {
                seqid: fields[0].to_string(),
                source: fields[1].to_string(),
                feature_type: fields[2].to_string(),
                start,
                end,
                score,
                strand: parse_strand(fields[6]),
                phase,
                attributes,
            }));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    #[test]
    fn test_gff_reader() {
        let data = "\
##gff-version 3
chr1\tEnsembl\tgene\t1001\t5000\t.\t+\t.\tID=gene1;Name=BRCA1
chr1\tEnsembl\texon\t1001\t1500\t.\t+\t0\tID=exon1;Parent=gene1
";
        let reader = BufReader::new(data.as_bytes());
        let records: Vec<_> = GffReader::new(reader).collect();
        assert_eq!(records.len(), 2);

        let gene = records[0].as_ref().unwrap();
        assert_eq!(gene.seqid, "chr1");
        assert_eq!(gene.feature_type, "gene");
        assert_eq!(gene.start, 1000); // 1-based 1001 -> 0-based 1000
        assert_eq!(gene.end, 5000);
    }

    #[test]
    fn test_gff_to_feature() {
        let data = "chr1\tTest\tgene\t1001\t5000\t.\t+\t.\tID=gene1\n";
        let reader = BufReader::new(data.as_bytes());
        let record = GffReader::new(reader).next().unwrap().unwrap();
        let feat = record.to_feature().unwrap();
        assert_eq!(feat.interval.start, 1000);
        assert_eq!(feat.id.as_deref(), Some("gene1"));
    }
}
