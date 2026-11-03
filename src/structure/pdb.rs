//! PDB and mmCIF 3D macromolecular structure parsers and coordinate primitives.

use crate::BioError;
use std::io::BufRead;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A 3D atom record in a macromolecular structure.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Atom3D {
    /// Atom serial number.
    pub serial: usize,
    /// Atom name (e.g. "CA", "CB", "N", "O").
    pub name: String,
    /// Residue name (e.g. "ALA", "GLY", "TRP").
    pub res_name: String,
    /// Chain identifier (e.g. 'A', 'B').
    pub chain_id: char,
    /// Residue sequence number.
    pub res_seq: i32,
    /// Cartesian X coordinate in Ångströms.
    pub x: f64,
    /// Cartesian Y coordinate in Ångströms.
    pub y: f64,
    /// Cartesian Z coordinate in Ångströms.
    pub z: f64,
    /// Crystallographic occupancy (typically 1.0).
    pub occupancy: f64,
    /// Temperature B-factor.
    pub b_factor: f64,
    /// Chemical element symbol (e.g. "C", "N", "O", "S").
    pub element: String,
}

impl Atom3D {
    /// Calculates Euclidean distance in Ångströms to another atom.
    pub fn distance(&self, other: &Atom3D) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Returns coordinates as `[x, y, z]`.
    pub fn coords(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
}

/// An amino acid or nucleotide residue comprising constituent atoms.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Residue3D {
    /// Three-letter residue code (e.g. "ALA").
    pub name: String,
    /// Residue sequence number.
    pub seq_num: i32,
    /// Constituent atoms.
    pub atoms: Vec<Atom3D>,
}

impl Residue3D {
    /// Returns the backbone Alpha Carbon (Cα) atom if present.
    pub fn ca_atom(&self) -> Option<&Atom3D> {
        self.atoms.iter().find(|a| a.name.trim() == "CA")
    }
}

/// A macromolecular chain (e.g. Chain A) composed of sequential residues.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Chain3D {
    /// Chain identifier character.
    pub id: char,
    /// Sequential residues in this chain.
    pub residues: Vec<Residue3D>,
}

/// A complete 3D protein or macromolecular structure.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ProteinStructure3D {
    /// Structure identifier (e.g. PDB ID "1CRN").
    pub id: String,
    /// Chains present in this structure.
    pub chains: Vec<Chain3D>,
}

impl ProteinStructure3D {
    /// Gathers all backbone Cα atoms across all chains.
    pub fn ca_atoms(&self) -> Vec<&Atom3D> {
        let mut cas = Vec::new();
        for chain in &self.chains {
            for res in &chain.residues {
                if let Some(ca) = res.ca_atom() {
                    cas.push(ca);
                }
            }
        }
        cas
    }
}

/// Parses standard Protein Data Bank (`.pdb`) format from a buffered stream.
pub fn parse_pdb<R: BufRead>(reader: R) -> crate::Result<ProteinStructure3D> {
    let mut structure_id = "UNKNOWN".to_string();
    let mut atoms: Vec<Atom3D> = Vec::new();

    for line_result in reader.lines() {
        let line = line_result.map_err(|e| BioError::ParseError(e.to_string()))?;
        if line.starts_with("HEADER") {
            if line.len() >= 66 {
                structure_id = line[62..66].trim().to_string();
            } else if let Some(last_tok) = line.split_whitespace().last() {
                if last_tok != "HEADER" {
                    structure_id = last_tok.to_string();
                }
            }
        } else if line.starts_with("ATOM") || line.starts_with("HETATM") {
            if line.len() < 54 {
                continue;
            }
            let serial: usize = line[6..11].trim().parse().unwrap_or(0);
            let name = line[12..16].trim().to_string();
            let res_name = line[17..20].trim().to_string();
            let chain_id = line[21..22].chars().next().unwrap_or('A');
            let res_seq: i32 = line[22..26].trim().parse().unwrap_or(0);

            let x: f64 = line[30..38].trim().parse().unwrap_or(0.0);
            let y: f64 = line[38..46].trim().parse().unwrap_or(0.0);
            let z: f64 = line[46..54].trim().parse().unwrap_or(0.0);

            let occupancy: f64 = if line.len() >= 60 {
                line[54..60].trim().parse().unwrap_or(1.0)
            } else {
                1.0
            };

            let b_factor: f64 = if line.len() >= 66 {
                line[60..66].trim().parse().unwrap_or(0.0)
            } else {
                0.0
            };

            let element = if line.len() >= 78 {
                line[76..78].trim().to_string()
            } else {
                name.chars().next().map(|c| c.to_string()).unwrap_or_else(|| "C".to_string())
            };

            atoms.push(Atom3D {
                serial,
                name,
                res_name,
                chain_id,
                res_seq,
                x,
                y,
                z,
                occupancy,
                b_factor,
                element,
            });
        }
    }

    // Group atoms into chains and residues
    let mut chains_map: std::collections::BTreeMap<char, std::collections::BTreeMap<i32, (String, Vec<Atom3D>)>> =
        std::collections::BTreeMap::new();

    for atom in atoms {
        chains_map
            .entry(atom.chain_id)
            .or_default()
            .entry(atom.res_seq)
            .or_insert_with(|| (atom.res_name.clone(), Vec::new()))
            .1
            .push(atom);
    }

    let mut chains = Vec::new();
    for (chain_id, res_map) in chains_map {
        let mut residues = Vec::new();
        for (seq_num, (name, atom_list)) in res_map {
            residues.push(Residue3D {
                name,
                seq_num,
                atoms: atom_list,
            });
        }
        chains.push(Chain3D {
            id: chain_id,
            residues,
        });
    }

    Ok(ProteinStructure3D {
        id: structure_id,
        chains,
    })
}

/// Parses mmCIF format records starting with `_atom_site` lines.
pub fn parse_mmcif<R: BufRead>(reader: R) -> crate::Result<ProteinStructure3D> {
    let mut atoms = Vec::new();
    let mut in_atom_site_loop = false;
    let mut col_map = std::collections::HashMap::new();
    let mut col_idx = 0;

    for line_result in reader.lines() {
        let line = line_result.map_err(|e| BioError::ParseError(e.to_string()))?;
        let trimmed = line.trim();

        if trimmed == "loop_" {
            in_atom_site_loop = false;
            col_map.clear();
            col_idx = 0;
            continue;
        }

        if trimmed.starts_with("_atom_site.") {
            in_atom_site_loop = true;
            let field = trimmed.strip_prefix("_atom_site.").unwrap_or("");
            col_map.insert(field.to_string(), col_idx);
            col_idx += 1;
            continue;
        }

        if in_atom_site_loop {
            if trimmed.starts_with('#') || trimmed.is_empty() {
                in_atom_site_loop = false;
                continue;
            }

            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() < col_map.len() {
                continue;
            }

            let group = parts.get(*col_map.get("group_PDB").unwrap_or(&0)).copied().unwrap_or("");
            if group != "ATOM" && group != "HETATM" {
                continue;
            }

            let id_idx = *col_map.get("id").unwrap_or(&1);
            let name_idx = *col_map.get("auth_atom_id").or_else(|| col_map.get("label_atom_id")).unwrap_or(&2);
            let res_idx = *col_map.get("auth_comp_id").or_else(|| col_map.get("label_comp_id")).unwrap_or(&3);
            let chain_idx = *col_map.get("auth_asym_id").or_else(|| col_map.get("label_asym_id")).unwrap_or(&4);
            let seq_idx = *col_map.get("auth_seq_id").or_else(|| col_map.get("label_seq_id")).unwrap_or(&5);
            let x_idx = *col_map.get("Cartn_x").unwrap_or(&6);
            let y_idx = *col_map.get("Cartn_y").unwrap_or(&7);
            let z_idx = *col_map.get("Cartn_z").unwrap_or(&8);
            let b_idx = *col_map.get("B_iso_or_equiv").unwrap_or(&9);

            let serial: usize = parts.get(id_idx).and_then(|s| s.parse().ok()).unwrap_or(0);
            let name = parts.get(name_idx).copied().unwrap_or("CA").to_string();
            let res_name = parts.get(res_idx).copied().unwrap_or("UNK").to_string();
            let chain_id = parts.get(chain_idx).and_then(|s| s.chars().next()).unwrap_or('A');
            let res_seq: i32 = parts.get(seq_idx).and_then(|s| s.parse().ok()).unwrap_or(0);
            let x: f64 = parts.get(x_idx).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let y: f64 = parts.get(y_idx).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let z: f64 = parts.get(z_idx).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let b_factor: f64 = parts.get(b_idx).and_then(|s| s.parse().ok()).unwrap_or(0.0);

            atoms.push(Atom3D {
                serial,
                name,
                res_name,
                chain_id,
                res_seq,
                x,
                y,
                z,
                occupancy: 1.0,
                b_factor,
                element: "C".to_string(),
            });
        }
    }

    // Build structure from atoms
    let mut chains_map: std::collections::BTreeMap<char, std::collections::BTreeMap<i32, (String, Vec<Atom3D>)>> =
        std::collections::BTreeMap::new();

    for atom in atoms {
        chains_map
            .entry(atom.chain_id)
            .or_default()
            .entry(atom.res_seq)
            .or_insert_with(|| (atom.res_name.clone(), Vec::new()))
            .1
            .push(atom);
    }

    let mut chains = Vec::new();
    for (chain_id, res_map) in chains_map {
        let mut residues = Vec::new();
        for (seq_num, (name, atom_list)) in res_map {
            residues.push(Residue3D {
                name,
                seq_num,
                atoms: atom_list,
            });
        }
        chains.push(Chain3D {
            id: chain_id,
            residues,
        });
    }

    Ok(ProteinStructure3D {
        id: "MMCIF".to_string(),
        chains,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_pdb_basic() {
        let pdb_content = "\
HEADER    SMALL PROTEIN                           1CRN\n\
ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N\n\
ATOM      2  CA  THR A   1      16.967  12.784   4.338  1.00 10.80           C\n\
ATOM      3  C   THR A   1      15.685  12.755   5.133  1.00  9.19           C\n\
ATOM      4  N   THR A   2      15.715  12.871   6.452  1.00  7.56           N\n\
ATOM      5  CA  THR A   2      14.542  12.868   7.309  1.00  7.02           C\n\
END\n";

        let reader = std::io::BufReader::new(pdb_content.as_bytes());
        let structure = parse_pdb(reader).unwrap();

        assert_eq!(structure.id, "1CRN");
        assert_eq!(structure.chains.len(), 1);
        let chain_a = &structure.chains[0];
        assert_eq!(chain_a.id, 'A');
        assert_eq!(chain_a.residues.len(), 2);

        let cas = structure.ca_atoms();
        assert_eq!(cas.len(), 2);
        assert_eq!(cas[0].res_seq, 1);
        assert_eq!(cas[1].res_seq, 2);

        let dist = cas[0].distance(cas[1]);
        assert!((dist - 3.8).abs() < 0.5); // Standard peptide Cα-Cα bond length ~3.8 Å
    }
}
