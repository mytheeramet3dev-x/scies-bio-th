use std::io::BufRead;

pub use crate::variant::{Variant, VcfRecord};

/// A streaming VCF parser (re-exported and mirrored from variant module for ergonomic IO access).
pub struct VcfReader<R> {
    inner: crate::variant::VcfReader<R>,
}

impl<R: BufRead> VcfReader<R> {
    pub fn new(reader: R) -> Self {
        Self { inner: crate::variant::VcfReader::new(reader) }
    }
}

impl<R: BufRead> Iterator for VcfReader<R> {
    type Item = crate::Result<VcfRecord>;
    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next()
    }
}
