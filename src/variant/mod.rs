use crate::genome::{Interval, Strand};
use crate::BioError;
use std::io::BufRead;

/// The classification of a genomic variant.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VariantType {
    Snp,
    Insertion,
    Deletion,
    Mnv,       // Multi-nucleotide variant
    Complex,
}

/// A genomic variant combining locus, reference, and alternate allele.
#[derive(Debug, Clone)]
pub struct Variant {
    pub interval: Interval,
    pub reference: Vec<u8>,
    pub alternate: Vec<u8>,
    pub variant_type: VariantType,
    pub quality: Option<f64>,
    pub id: Option<String>,
}

impl Variant {
    pub fn new(
        chrom: impl Into<String>,
        pos: u64,          // 0-based
        reference: Vec<u8>,
        alternate: Vec<u8>,
    ) -> crate::Result<Self> {
        let end = pos + reference.len() as u64;
        let interval = Interval::new(chrom, pos, end, Strand::Unknown)?;

        let variant_type = if reference.len() == 1 && alternate.len() == 1 {
            VariantType::Snp
        } else if reference.len() == 1 {
            VariantType::Insertion
        } else if alternate.len() == 1 {
            VariantType::Deletion
        } else if reference.len() == alternate.len() {
            VariantType::Mnv
        } else {
            VariantType::Complex
        };

        Ok(Self {
            interval,
            reference,
            alternate,
            variant_type,
            quality: None,
            id: None,
        })
    }
}

/// A simplified VCF record.
#[derive(Debug, Clone)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: u64,       // 0-based
    pub id: Option<String>,
    pub reference: Vec<u8>,
    pub alternate: Vec<u8>,
    pub quality: Option<f64>,
    pub filter: String,
    pub info: String,
}

impl VcfRecord {
    /// Convert a VCF record into a `Variant`.
    pub fn to_variant(&self) -> crate::Result<Variant> {
        let mut v = Variant::new(
            &self.chrom,
            self.pos,
            self.reference.clone(),
            self.alternate.clone(),
        )?;
        v.quality = self.quality;
        v.id = self.id.clone();
        Ok(v)
    }
}

/// A streaming VCF parser.
pub struct VcfReader<R> {
    reader: R,
}

impl<R: BufRead> VcfReader<R> {
    pub fn new(reader: R) -> Self {
        Self { reader }
    }
}

impl<R: BufRead> Iterator for VcfReader<R> {
    type Item = crate::Result<VcfRecord>;

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
                continue;
            }

            let fields: Vec<&str> = trimmed.splitn(9, '\t').collect();
            if fields.len() < 5 {
                return Some(Err(BioError::ParseError(
                    format!("VCF: expected at least 5 fields, got {}", fields.len())
                )));
            }

            // VCF POS is 1-based; convert to 0-based
            let pos = match fields[1].parse::<u64>() {
                Ok(v) => v.saturating_sub(1),
                Err(_) => return Some(Err(BioError::ParseError("Invalid VCF POS".into()))),
            };

            let id = match fields[2] {
                "." => None,
                s => Some(s.to_string()),
            };

            let quality = if fields.len() > 5 {
                match fields[5] {
                    "." => None,
                    s => s.parse::<f64>().ok(),
                }
            } else {
                None
            };

            let filter = if fields.len() > 6 { fields[6].to_string() } else { ".".to_string() };
            let info = if fields.len() > 7 { fields[7].to_string() } else { ".".to_string() };

            return Some(Ok(VcfRecord {
                chrom: fields[0].to_string(),
                pos,
                id,
                reference: fields[3].as_bytes().to_vec(),
                alternate: fields[4].as_bytes().to_vec(),
                quality,
                filter,
                info,
            }));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    #[test]
    fn test_vcf_reader() {
        let data = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t12345\t.\tA\tG\t50.0\tPASS\t.
chr1\t67890\trs123\tACT\tA\t30.0\tPASS\t.
";
        let reader = BufReader::new(data.as_bytes());
        let records: Vec<_> = VcfReader::new(reader).collect();
        assert_eq!(records.len(), 2);

        let snp = records[0].as_ref().unwrap();
        assert_eq!(snp.pos, 12344); // 0-based
        assert_eq!(snp.reference, b"A");
        assert_eq!(snp.alternate, b"G");

        let del = records[1].as_ref().unwrap();
        assert_eq!(del.id.as_deref(), Some("rs123"));
    }

    #[test]
    fn test_variant_type() {
        let snp = Variant::new("chr1", 100, b"A".to_vec(), b"G".to_vec()).unwrap();
        assert_eq!(snp.variant_type, VariantType::Snp);

        let ins = Variant::new("chr1", 200, b"A".to_vec(), b"AGT".to_vec()).unwrap();
        assert_eq!(ins.variant_type, VariantType::Insertion);

        let del = Variant::new("chr1", 300, b"AGT".to_vec(), b"A".to_vec()).unwrap();
        assert_eq!(del.variant_type, VariantType::Deletion);
    }
}
