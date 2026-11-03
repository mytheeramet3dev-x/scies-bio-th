use std::cmp::max;

pub struct AlignmentResult {
    pub score: i32,
    pub alignment1: String,
    pub alignment2: String,
}

pub struct NeedlemanWunsch {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_penalty: i32,
}

impl NeedlemanWunsch {
    pub fn new(match_score: i32, mismatch_score: i32, gap_penalty: i32) -> Self {
        Self { match_score, mismatch_score, gap_penalty }
    }

    pub fn align(&self, seq1: &[u8], seq2: &[u8]) -> AlignmentResult {
        let n = seq1.len();
        let m = seq2.len();
        let cols = m + 1;
        let mut matrix = vec![0i32; (n + 1) * cols];

        // Initialization
        for (i, row) in matrix.chunks_exact_mut(cols).enumerate().take(n + 1) {
            row[0] = i as i32 * self.gap_penalty;
        }
        for (j, cell) in matrix[..cols].iter_mut().enumerate() {
            *cell = j as i32 * self.gap_penalty;
        }

        // DP filling
        for i in 1..=n {
            let row_offset = i * cols;
            let prev_row_offset = (i - 1) * cols;
            for j in 1..=m {
                let match_mismatch = if seq1[i - 1] == seq2[j - 1] {
                    self.match_score
                } else {
                    self.mismatch_score
                };

                let diag = matrix[prev_row_offset + j - 1] + match_mismatch;
                let up = matrix[prev_row_offset + j] + self.gap_penalty;
                let left = matrix[row_offset + j - 1] + self.gap_penalty;

                matrix[row_offset + j] = diag.max(up).max(left);
            }
        }

        // Traceback
        let mut align1 = Vec::new();
        let mut align2 = Vec::new();
        let mut i = n;
        let mut j = m;
        let score = matrix[n * cols + m];

        while i > 0 || j > 0 {
            let curr_offset = i * cols + j;
            if i > 0 && j > 0 {
                let current_score = matrix[curr_offset];
                let match_mismatch = if seq1[i - 1] == seq2[j - 1] {
                    self.match_score
                } else {
                    self.mismatch_score
                };

                if current_score == matrix[(i - 1) * cols + (j - 1)] + match_mismatch {
                    align1.push(seq1[i - 1]);
                    align2.push(seq2[j - 1]);
                    i -= 1;
                    j -= 1;
                    continue;
                }
            }
            if i > 0 && (j == 0 || matrix[curr_offset] == matrix[(i - 1) * cols + j] + self.gap_penalty) {
                align1.push(seq1[i - 1]);
                align2.push(b'-');
                i -= 1;
            } else {
                align1.push(b'-');
                align2.push(seq2[j - 1]);
                j -= 1;
            }
        }

        align1.reverse();
        align2.reverse();

        AlignmentResult {
            score,
            alignment1: String::from_utf8_lossy(&align1).into_owned(),
            alignment2: String::from_utf8_lossy(&align2).into_owned(),
        }
    }
}

pub struct SmithWaterman {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_penalty: i32,
}

impl SmithWaterman {
    pub fn new(match_score: i32, mismatch_score: i32, gap_penalty: i32) -> Self {
        Self { match_score, mismatch_score, gap_penalty }
    }

    pub fn align(&self, seq1: &[u8], seq2: &[u8]) -> AlignmentResult {
        let n = seq1.len();
        let m = seq2.len();
        let cols = m + 1;
        let mut matrix = vec![0i32; (n + 1) * cols];
        let mut max_score = 0;
        let mut max_i = 0;
        let mut max_j = 0;

        // DP filling
        for i in 1..=n {
            let row_offset = i * cols;
            let prev_row_offset = (i - 1) * cols;
            for j in 1..=m {
                let match_mismatch = if seq1[i - 1] == seq2[j - 1] {
                    self.match_score
                } else {
                    self.mismatch_score
                };

                let diag = matrix[prev_row_offset + j - 1] + match_mismatch;
                let up = matrix[prev_row_offset + j] + self.gap_penalty;
                let left = matrix[row_offset + j - 1] + self.gap_penalty;

                let score = max(0, diag.max(up).max(left));
                matrix[row_offset + j] = score;

                if score > max_score {
                    max_score = score;
                    max_i = i;
                    max_j = j;
                }
            }
        }

        // Traceback from max_i, max_j
        let mut align1 = Vec::new();
        let mut align2 = Vec::new();
        let mut i = max_i;
        let mut j = max_j;

        while i > 0 && j > 0 && matrix[i * cols + j] > 0 {
            let current_score = matrix[i * cols + j];
            let match_mismatch = if seq1[i - 1] == seq2[j - 1] {
                self.match_score
            } else {
                self.mismatch_score
            };

            if current_score == matrix[(i - 1) * cols + (j - 1)] + match_mismatch {
                align1.push(seq1[i - 1]);
                align2.push(seq2[j - 1]);
                i -= 1;
                j -= 1;
            } else if current_score == matrix[(i - 1) * cols + j] + self.gap_penalty {
                align1.push(seq1[i - 1]);
                align2.push(b'-');
                i -= 1;
            } else {
                align1.push(b'-');
                align2.push(seq2[j - 1]);
                j -= 1;
            }
        }

        align1.reverse();
        align2.reverse();

        AlignmentResult {
            score: max_score,
            alignment1: String::from_utf8_lossy(&align1).into_owned(),
            alignment2: String::from_utf8_lossy(&align2).into_owned(),
        }
    }
}

pub fn levenshtein_distance(seq1: &[u8], seq2: &[u8]) -> usize {
    let n = seq1.len();
    let m = seq2.len();
    if n == 0 { return m; }
    if m == 0 { return n; }

    let mut dp = vec![0; m + 1];
    for (j, cell) in dp.iter_mut().enumerate() {
        *cell = j;
    }

    for i in 1..=n {
        let mut prev_diag = dp[0];
        dp[0] = i;
        for j in 1..=m {
            let temp = dp[j];
            if seq1[i - 1] == seq2[j - 1] {
                dp[j] = prev_diag;
            } else {
                dp[j] = 1 + dp[j].min(dp[j - 1]).min(prev_diag);
            }
            prev_diag = temp;
        }
    }
    dp[m]
}

pub fn hamming_distance(seq1: &[u8], seq2: &[u8]) -> crate::Result<usize> {
    if seq1.len() != seq2.len() {
        return Err(crate::BioError::AlignmentMismatch);
    }
    let dist = seq1.iter().zip(seq2.iter()).filter(|(a, b)| a != b).count();
    Ok(dist)
}

const NEG_INF: i32 = i32::MIN / 4;

/// Banded Needleman-Wunsch global alignment with band radius `k`.
///
/// Only evaluates dynamic programming matrix cells where $|i - j| \le k$.
/// Complexity is reduced from $O(M \cdot N)$ to $O(k \cdot \min(M, N))$.
pub fn banded_needleman_wunsch(
    seq1: &[u8],
    seq2: &[u8],
    match_score: i32,
    mismatch_score: i32,
    gap_penalty: i32,
    k: usize,
) -> AlignmentResult {
    let n = seq1.len();
    let m = seq2.len();
    let cols = m + 1;
    let mut matrix = vec![NEG_INF; (n + 1) * cols];

    matrix[0] = 0;
    for i in 1..=(n.min(k)) {
        matrix[i * cols] = i as i32 * gap_penalty;
    }
    for (j, cell) in matrix[..=(m.min(k))].iter_mut().enumerate().skip(1) {
        *cell = j as i32 * gap_penalty;
    }

    let k_signed = k as isize;

    for i in 1..=n {
        let i_s = i as isize;
        let j_start = (i_s - k_signed).max(1) as usize;
        let j_end = (i_s + k_signed).min(m as isize) as usize;

        let row_offset = i * cols;
        let prev_row_offset = (i - 1) * cols;

        for j in j_start..=j_end {
            let match_mismatch = if seq1[i - 1] == seq2[j - 1] {
                match_score
            } else {
                mismatch_score
            };

            let diag = matrix[prev_row_offset + j - 1] + match_mismatch;
            let up = matrix[prev_row_offset + j] + gap_penalty;
            let left = matrix[row_offset + j - 1] + gap_penalty;

            matrix[row_offset + j] = diag.max(up).max(left);
        }
    }

    // Traceback
    let mut align1 = Vec::new();
    let mut align2 = Vec::new();
    let mut i = n;
    let mut j = m;
    let score = matrix[n * cols + m];

    while i > 0 || j > 0 {
        let curr_offset = i * cols + j;
        if i > 0 && j > 0 {
            let current_score = matrix[curr_offset];
            let match_mismatch = if seq1[i - 1] == seq2[j - 1] {
                match_score
            } else {
                mismatch_score
            };

            if current_score == matrix[(i - 1) * cols + (j - 1)] + match_mismatch {
                align1.push(seq1[i - 1]);
                align2.push(seq2[j - 1]);
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if i > 0 && (j == 0 || matrix[curr_offset] == matrix[(i - 1) * cols + j] + gap_penalty) {
            align1.push(seq1[i - 1]);
            align2.push(b'-');
            i -= 1;
        } else if j > 0 {
            align1.push(b'-');
            align2.push(seq2[j - 1]);
            j -= 1;
        } else {
            break;
        }
    }

    align1.reverse();
    align2.reverse();

    AlignmentResult {
        score,
        alignment1: String::from_utf8_lossy(&align1).into_owned(),
        alignment2: String::from_utf8_lossy(&align2).into_owned(),
    }
}

/// Banded Smith-Waterman local alignment with band radius `k`.
pub fn banded_smith_waterman(
    seq1: &[u8],
    seq2: &[u8],
    match_score: i32,
    mismatch_score: i32,
    gap_penalty: i32,
    k: usize,
) -> AlignmentResult {
    let n = seq1.len();
    let m = seq2.len();
    let cols = m + 1;
    let mut matrix = vec![0i32; (n + 1) * cols];

    let mut max_score = 0;
    let mut max_i = 0;
    let mut max_j = 0;

    let k_signed = k as isize;

    for i in 1..=n {
        let i_s = i as isize;
        let j_start = (i_s - k_signed).max(1) as usize;
        let j_end = (i_s + k_signed).min(m as isize) as usize;

        let row_offset = i * cols;
        let prev_row_offset = (i - 1) * cols;

        for j in j_start..=j_end {
            let match_mismatch = if seq1[i - 1] == seq2[j - 1] {
                match_score
            } else {
                mismatch_score
            };

            let diag = matrix[prev_row_offset + j - 1] + match_mismatch;
            let up = matrix[prev_row_offset + j] + gap_penalty;
            let left = matrix[row_offset + j - 1] + gap_penalty;

            let score = 0.max(diag.max(up).max(left));
            matrix[row_offset + j] = score;

            if score > max_score {
                max_score = score;
                max_i = i;
                max_j = j;
            }
        }
    }

    let mut align1 = Vec::new();
    let mut align2 = Vec::new();
    let mut i = max_i;
    let mut j = max_j;

    while i > 0 && j > 0 && matrix[i * cols + j] > 0 {
        let current_score = matrix[i * cols + j];
        let match_mismatch = if seq1[i - 1] == seq2[j - 1] {
            match_score
        } else {
            mismatch_score
        };

        if current_score == matrix[(i - 1) * cols + (j - 1)] + match_mismatch {
            align1.push(seq1[i - 1]);
            align2.push(seq2[j - 1]);
            i -= 1;
            j -= 1;
        } else if current_score == matrix[(i - 1) * cols + j] + gap_penalty {
            align1.push(seq1[i - 1]);
            align2.push(b'-');
            i -= 1;
        } else {
            align1.push(b'-');
            align2.push(seq2[j - 1]);
            j -= 1;
        }
    }

    align1.reverse();
    align2.reverse();

    AlignmentResult {
        score: max_score,
        alignment1: String::from_utf8_lossy(&align1).into_owned(),
        alignment2: String::from_utf8_lossy(&align2).into_owned(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_needleman_wunsch() {
        let nw = NeedlemanWunsch::new(1, -1, -1);
        let res = nw.align(b"GCATGCU", b"GATTACA");
        assert_eq!(res.score, 0); 
    }

    #[test]
    fn test_smith_waterman() {
        let sw = SmithWaterman::new(2, -1, -1);
        let res = sw.align(b"TGTTACGG", b"GGTTGACTA");
        assert!(res.score > 0);
    }

    #[test]
    fn test_levenshtein() {
        assert_eq!(levenshtein_distance(b"kitten", b"sitting"), 3);
        assert_eq!(levenshtein_distance(b"flaw", b"lawn"), 2);
    }

    #[test]
    fn test_hamming() {
        assert_eq!(hamming_distance(b"karolin", b"kathrin").unwrap(), 3);
        assert!(hamming_distance(b"karolin", b"kathrinn").is_err());
    }

    #[test]
    fn test_banded_needleman_wunsch() {
        let seq1 = b"GCATGCU";
        let seq2 = b"GATTACA";
        let res_banded = banded_needleman_wunsch(seq1, seq2, 1, -1, -1, 3);
        let nw = NeedlemanWunsch::new(1, -1, -1);
        let res_full = nw.align(seq1, seq2);
        assert_eq!(res_banded.score, res_full.score);
    }

    #[test]
    fn test_banded_smith_waterman() {
        let seq1 = b"TGTTACGG";
        let seq2 = b"GGTTGACTA";
        let res_banded = banded_smith_waterman(seq1, seq2, 2, -1, -1, 4);
        let sw = SmithWaterman::new(2, -1, -1);
        let res_full = sw.align(seq1, seq2);
        assert_eq!(res_banded.score, res_full.score);
    }
}
