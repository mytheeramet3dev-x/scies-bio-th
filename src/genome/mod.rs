use std::cmp::Ordering;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A strand of a genomic feature.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

/// A genomic interval [start, end) on a given chromosome/contig.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Interval {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
}

impl Interval {
    pub fn new(chrom: impl Into<String>, start: u64, end: u64, strand: Strand) -> crate::Result<Self> {
        let chrom = chrom.into();
        if start >= end {
            return Err(crate::BioError::ParseError(format!(
                "Interval start ({}) must be < end ({})",
                start, end
            )));
        }
        Ok(Self { chrom, start, end, strand })
    }

    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    pub fn is_empty(&self) -> bool {
        self.start >= self.end
    }

    pub fn overlaps(&self, other: &Interval) -> bool {
        self.chrom == other.chrom && self.start < other.end && other.start < self.end
    }

    pub fn contains(&self, pos: u64) -> bool {
        pos >= self.start && pos < self.end
    }
}

impl PartialOrd for Interval {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Interval {
    fn cmp(&self, other: &Self) -> Ordering {
        self.chrom
            .cmp(&other.chrom)
            .then(self.start.cmp(&other.start))
            .then(self.end.cmp(&other.end))
    }
}

/// A simple sorted-list interval store that supports overlap queries.
/// For v0.2 we use a sorted Vec with binary search; a full augmented tree
/// can be introduced in v0.3 when performance profiling demands it.
pub struct IntervalTree<T> {
    entries: Vec<(Interval, T)>,
    sorted: bool,
}

impl<T: Clone> IntervalTree<T> {
    pub fn new() -> Self {
        Self { entries: Vec::new(), sorted: false }
    }

    pub fn insert(&mut self, interval: Interval, data: T) {
        self.entries.push((interval, data));
        self.sorted = false;
    }

    pub fn sort(&mut self) {
        self.entries.sort_by(|(a, _), (b, _)| a.cmp(b));
        self.sorted = true;
    }

    pub fn query_overlapping(&mut self, query: &Interval) -> Vec<&T> {
        if !self.sorted {
            self.sort();
        }
        self.entries
            .iter()
            .filter(|(iv, _)| iv.overlaps(query))
            .map(|(_, data)| data)
            .collect()
    }
}

impl<T: Clone> Default for IntervalTree<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interval_overlap() {
        let a = Interval::new("chr1", 100, 200, Strand::Forward).unwrap();
        let b = Interval::new("chr1", 150, 250, Strand::Forward).unwrap();
        let c = Interval::new("chr1", 300, 400, Strand::Forward).unwrap();
        assert!(a.overlaps(&b));
        assert!(!a.overlaps(&c));
    }

    #[test]
    fn test_interval_contains() {
        let iv = Interval::new("chr1", 100, 200, Strand::Forward).unwrap();
        assert!(iv.contains(100));
        assert!(iv.contains(150));
        assert!(!iv.contains(200));
    }

    #[test]
    fn test_interval_invalid() {
        assert!(Interval::new("chr1", 200, 100, Strand::Forward).is_err());
    }

    #[test]
    fn test_interval_tree_query() {
        let mut tree: IntervalTree<&str> = IntervalTree::new();
        tree.insert(Interval::new("chr1", 100, 200, Strand::Forward).unwrap(), "geneA");
        tree.insert(Interval::new("chr1", 300, 500, Strand::Forward).unwrap(), "geneB");

        let q = Interval::new("chr1", 150, 350, Strand::Unknown).unwrap();
        let hits = tree.query_overlapping(&q);
        assert_eq!(hits.len(), 2);
    }
}
