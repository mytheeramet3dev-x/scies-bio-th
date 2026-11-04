use crate::BioError;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

pub mod sha256;
pub use sha256::{compute_sha256_hex, Sha256};

/// A DNA sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Dna {
    pub sequence: Vec<u8>,
}

impl Dna {
    /// Creates a new `Dna` sequence.
    /// Validates that all characters are valid DNA bases (A, C, G, T, N).
    /// Both uppercase and lowercase are accepted, but they are converted to uppercase.
    pub fn new(seq: &[u8]) -> crate::Result<Self> {
        let mut sequence = Vec::with_capacity(seq.len());
        for &byte in seq {
            let b = byte.to_ascii_uppercase();
            match b {
                b'A' | b'C' | b'G' | b'T' | b'N' => sequence.push(b),
                _ => return Err(BioError::InvalidCharacter(byte as char)),
            }
        }
        Ok(Self { sequence })
    }

    /// Returns the reverse complement of the DNA sequence.
    pub fn reverse_complement(&self) -> Self {
        let mut rc = Vec::with_capacity(self.sequence.len());
        for &byte in self.sequence.iter().rev() {
            let comp = match byte {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                b'N' => b'N',
                _ => unreachable!(),
            };
            rc.push(comp);
        }
        Self { sequence: rc }
    }

    /// Calculates the GC content of the DNA sequence as a percentage (0.0 to 100.0).
    pub fn gc_content(&self) -> f64 {
        if self.sequence.is_empty() {
            return 0.0;
        }
        let gc_count = self.sequence.iter().filter(|&&b| b == b'G' || b == b'C').count();
        (gc_count as f64 / self.sequence.len() as f64) * 100.0
    }

    /// Returns the underlying byte slice.
    pub fn as_bytes(&self) -> &[u8] {
        &self.sequence
    }
}

/// An RNA sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Rna {
    pub sequence: Vec<u8>,
}

impl Rna {
    /// Creates a new `Rna` sequence.
    /// Validates that all characters are valid RNA bases (A, C, G, U, N).
    pub fn new(seq: &[u8]) -> crate::Result<Self> {
        let mut sequence = Vec::with_capacity(seq.len());
        for &byte in seq {
            let b = byte.to_ascii_uppercase();
            match b {
                b'A' | b'C' | b'G' | b'U' | b'N' => sequence.push(b),
                _ => return Err(BioError::InvalidCharacter(byte as char)),
            }
        }
        Ok(Self { sequence })
    }

    /// Calculates the GC content of the RNA sequence.
    pub fn gc_content(&self) -> f64 {
        if self.sequence.is_empty() {
            return 0.0;
        }
        let gc_count = self.sequence.iter().filter(|&&b| b == b'G' || b == b'C').count();
        gc_count as f64 / self.sequence.len() as f64
    }

    /// Returns the underlying byte slice.
    pub fn as_bytes(&self) -> &[u8] {
        &self.sequence
    }
}

/// A Protein sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Protein {
    pub sequence: Vec<u8>,
}

impl Protein {
    /// Creates a new `Protein` sequence.
    /// Validates standard amino acid characters.
    pub fn new(seq: &[u8]) -> crate::Result<Self> {
        let mut sequence = Vec::with_capacity(seq.len());
        for &byte in seq {
            let b = byte.to_ascii_uppercase();
            match b {
                b'A' | b'R' | b'N' | b'D' | b'C' | b'E' | b'Q' | b'G' | b'H' | b'I' |
                b'L' | b'K' | b'M' | b'F' | b'P' | b'S' | b'T' | b'W' | b'Y' | b'V' | b'*' | b'X' => sequence.push(b),
                _ => return Err(BioError::InvalidCharacter(byte as char)),
            }
        }
        Ok(Self { sequence })
    }

    /// Calculates the GC content of the RNA sequence.
    pub fn gc_content(&self) -> f64 {
        if self.sequence.is_empty() {
            return 0.0;
        }
        let gc_count = self.sequence.iter().filter(|&&b| b == b'G' || b == b'C').count();
        gc_count as f64 / self.sequence.len() as f64
    }

    /// Returns the underlying byte slice.
    pub fn as_bytes(&self) -> &[u8] {
        &self.sequence
    }
}

/// A fixed-size K-mer sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Kmer<const N: usize> {
    pub sequence: [u8; N],
}

impl<const N: usize> Kmer<N> {
    /// Creates a new `Kmer`.
    pub fn new(seq: &[u8; N]) -> Self {
        let mut sequence = [0u8; N];
        sequence.copy_from_slice(seq);
        Self { sequence }
    }
}

#[cfg(feature = "serde")]
impl<const N: usize> serde::Serialize for Kmer<N> {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_bytes(&self.sequence)
    }
}

#[cfg(feature = "serde")]
impl<'de, const N: usize> serde::Deserialize<'de> for Kmer<N> {
    fn deserialize<D>(deserializer: D) -> std::result::Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        struct KmerVisitor<const N: usize>;
        impl<'de, const N: usize> serde::de::Visitor<'de> for KmerVisitor<N> {
            type Value = Kmer<N>;

            fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(formatter, "a byte array of length {}", N)
            }

            fn visit_bytes<E>(self, v: &[u8]) -> std::result::Result<Kmer<N>, E>
            where
                E: serde::de::Error,
            {
                if v.len() != N {
                    return Err(E::custom(format!("expected length {}, got {}", N, v.len())));
                }
                let mut arr = [0u8; N];
                arr.copy_from_slice(v);
                Ok(Kmer { sequence: arr })
            }

            fn visit_seq<A>(self, mut seq: A) -> std::result::Result<Kmer<N>, A::Error>
            where
                A: serde::de::SeqAccess<'de>,
            {
                let mut arr = [0u8; N];
                for item in arr.iter_mut().take(N) {
                    *item = seq
                        .next_element()?
                        .ok_or_else(|| serde::de::Error::custom("insufficient elements"))?;
                }
                Ok(Kmer { sequence: arr })
            }
        }
        deserializer.deserialize_bytes(KmerVisitor::<N>)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dna_new_valid() {
        let dna = Dna::new(b"ACGTNacgtn").unwrap();
        assert_eq!(dna.sequence, b"ACGTNACGTN");
    }

    #[test]
    fn test_dna_new_invalid() {
        let res = Dna::new(b"ACGTU");
        assert!(res.is_err());
    }

    #[test]
    fn test_reverse_complement() {
        let dna = Dna::new(b"ATCGN").unwrap();
        let rc = dna.reverse_complement();
        assert_eq!(rc.sequence, b"NCGAT");
    }

    #[test]
    fn test_gc_content() {
        let dna = Dna::new(b"ATGCATGC").unwrap(); // 4 GC out of 8 = 50.0%
        assert_eq!(dna.gc_content(), 50.0);
    }
    
    #[test]
    fn test_rna_new() {
        let rna = Rna::new(b"ACGU").unwrap();
        assert_eq!(rna.sequence, b"ACGU");
        assert!(Rna::new(b"T").is_err());
    }

    #[test]
    fn test_protein_new() {
        let protein = Protein::new(b"ARNDCEQGHILKMFPSTWYVX*").unwrap();
        assert_eq!(protein.sequence.len(), 22);
        assert!(Protein::new(b"J").is_err());
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_kmer_serde() {
        let kmer: Kmer<4> = Kmer::new(b"ATGC");
        let serialized = serde_json::to_string(&kmer).unwrap();
        let deserialized: Kmer<4> = serde_json::from_str(&serialized).unwrap();
        assert_eq!(kmer, deserialized);
    }
}
