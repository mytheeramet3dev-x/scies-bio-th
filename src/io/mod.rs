pub mod gff;
pub mod vcf;

use crate::BioError;
use std::io::BufRead;

pub struct FastaRecord {
    pub id: String,
    pub sequence: Vec<u8>,
}

pub struct FastaReader<R> {
    reader: R,
    next_id: Option<String>,
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(mut reader: R) -> crate::Result<Self> {
        let mut first_line = String::new();
        let bytes_read = reader.read_line(&mut first_line)
            .map_err(|e| BioError::ParseError(e.to_string()))?;
        
        if bytes_read == 0 {
            return Ok(Self { reader, next_id: None });
        }
        
        let trimmed = first_line.trim_end();
        let next_id = match trimmed.strip_prefix('>') {
            Some(id) => Some(id.to_string()),
            None => return Err(BioError::ParseError("FASTA must start with '>'".to_string())),
        };
        
        Ok(Self {
            reader,
            next_id,
        })
    }
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = crate::Result<FastaRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let id = self.next_id.take()?;
        let mut sequence = Vec::new();
        let mut next_id = None;
        let mut line = String::new();

        loop {
            line.clear();
            match self.reader.read_line(&mut line) {
                Ok(0) => break, // EOF
                Ok(_) => {
                    let trimmed = line.trim_end();
                    if let Some(stripped) = trimmed.strip_prefix('>') {
                        next_id = Some(stripped.to_string());
                        break;
                    } else if !trimmed.is_empty() {
                        sequence.extend_from_slice(trimmed.as_bytes());
                    }
                }
                Err(e) => return Some(Err(BioError::ParseError(e.to_string()))),
            }
        }

        self.next_id = next_id;
        Some(Ok(FastaRecord { id, sequence }))
    }
}

pub struct FastqRecord {
    pub id: String,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
}

pub struct FastqReader<R> {
    reader: R,
}

impl<R: BufRead> FastqReader<R> {
    pub fn new(reader: R) -> Self {
        Self { reader }
    }
}

impl<R: BufRead> Iterator for FastqReader<R> {
    type Item = crate::Result<FastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut id_line = String::new();
        let mut seq_line = String::new();
        let mut plus_line = String::new();
        let mut qual_line = String::new();

        match self.reader.read_line(&mut id_line) {
            Ok(0) => return None, // EOF
            Err(e) => return Some(Err(BioError::ParseError(e.to_string()))),
            _ => {}
        }

        let id = id_line.trim_end();
        if !id.starts_with('@') {
            return Some(Err(BioError::ParseError("FASTQ record must start with '@'".to_string())));
        }
        let id = id[1..].to_string();

        if let Err(e) = self.reader.read_line(&mut seq_line) {
            return Some(Err(BioError::ParseError(e.to_string())));
        }
        if let Err(e) = self.reader.read_line(&mut plus_line) {
            return Some(Err(BioError::ParseError(e.to_string())));
        }
        if let Err(e) = self.reader.read_line(&mut qual_line) {
            return Some(Err(BioError::ParseError(e.to_string())));
        }
        
        let seq = seq_line.trim_end().as_bytes().to_vec();
        let plus = plus_line.trim_end();
        let qual = qual_line.trim_end().as_bytes().to_vec();

        if !plus.starts_with('+') {
            return Some(Err(BioError::ParseError("FASTQ separator must start with '+'".to_string())));
        }
        if seq.len() != qual.len() {
            return Some(Err(BioError::ParseError("FASTQ sequence and quality lengths differ".to_string())));
        }

        Some(Ok(FastqRecord {
            id,
            sequence: seq,
            quality: qual,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    #[test]
    fn test_fasta_reader() {
        let data = ">seq1\nATGC\n>seq2\nCGTA\n";
        let reader = BufReader::new(data.as_bytes());
        let mut fasta = FastaReader::new(reader).unwrap();
        
        let r1 = fasta.next().unwrap().unwrap();
        assert_eq!(r1.id, "seq1");
        assert_eq!(r1.sequence, b"ATGC");

        let r2 = fasta.next().unwrap().unwrap();
        assert_eq!(r2.id, "seq2");
        assert_eq!(r2.sequence, b"CGTA");

        assert!(fasta.next().is_none());
    }

    #[test]
    fn test_fasta_multiline() {
        let data = ">seq1\nAT\nGC\n>seq2\nC\nG\nTA";
        let reader = BufReader::new(data.as_bytes());
        let mut fasta = FastaReader::new(reader).unwrap();
        
        let r1 = fasta.next().unwrap().unwrap();
        assert_eq!(r1.id, "seq1");
        assert_eq!(r1.sequence, b"ATGC");

        let r2 = fasta.next().unwrap().unwrap();
        assert_eq!(r2.sequence, b"CGTA");
    }

    #[test]
    fn test_fastq_reader() {
        let data = "@seq1\nATGC\n+\n!!!!\n@seq2\nCGTA\n+\n~~~~";
        let reader = BufReader::new(data.as_bytes());
        let mut fastq = FastqReader::new(reader);

        let r1 = fastq.next().unwrap().unwrap();
        assert_eq!(r1.id, "seq1");
        assert_eq!(r1.sequence, b"ATGC");
        assert_eq!(r1.quality, b"!!!!");

        let r2 = fastq.next().unwrap().unwrap();
        assert_eq!(r2.id, "seq2");
        assert_eq!(r2.sequence, b"CGTA");
        assert_eq!(r2.quality, b"~~~~");

        assert!(fastq.next().is_none());
    }
}
