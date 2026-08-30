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
        let mut score_matrix = vec![vec![0; m + 1]; n + 1];

        // Initialization
        for i in 0..=n {
            score_matrix[i][0] = i as i32 * self.gap_penalty;
        }
        for j in 0..=m {
            score_matrix[0][j] = j as i32 * self.gap_penalty;
        }

        // DP filling
        for i in 1..=n {
            for j in 1..=m {
                let match_mismatch = if seq1[i - 1] == seq2[j - 1] {
                    self.match_score
                } else {
                    self.mismatch_score
                };

                let diag = score_matrix[i - 1][j - 1] + match_mismatch;
                let up = score_matrix[i - 1][j] + self.gap_penalty;
                let left = score_matrix[i][j - 1] + self.gap_penalty;

                score_matrix[i][j] = diag.max(up).max(left);
            }
        }

        // Traceback
        let mut align1 = Vec::new();
        let mut align2 = Vec::new();
        let mut i = n;
        let mut j = m;
        let score = score_matrix[n][m];

        while i > 0 || j > 0 {
            if i > 0 && j > 0 {
                let current_score = score_matrix[i][j];
                let match_mismatch = if seq1[i - 1] == seq2[j - 1] {
                    self.match_score
                } else {
                    self.mismatch_score
                };

                if current_score == score_matrix[i - 1][j - 1] + match_mismatch {
                    align1.push(seq1[i - 1]);
                    align2.push(seq2[j - 1]);
                    i -= 1;
                    j -= 1;
                    continue;
                }
            }
            if i > 0 && (j == 0 || score_matrix[i][j] == score_matrix[i - 1][j] + self.gap_penalty) {
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
        let mut score_matrix = vec![vec![0; m + 1]; n + 1];
        let mut max_score = 0;
        let mut max_i = 0;
        let mut max_j = 0;

        // DP filling (initialization is 0, which is handled by vec![vec![0; ...]])
        for i in 1..=n {
            for j in 1..=m {
                let match_mismatch = if seq1[i - 1] == seq2[j - 1] {
                    self.match_score
                } else {
                    self.mismatch_score
                };

                let diag = score_matrix[i - 1][j - 1] + match_mismatch;
                let up = score_matrix[i - 1][j] + self.gap_penalty;
                let left = score_matrix[i][j - 1] + self.gap_penalty;

                let score = max(0, diag.max(up).max(left));
                score_matrix[i][j] = score;

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

        while i > 0 && j > 0 && score_matrix[i][j] > 0 {
            let current_score = score_matrix[i][j];
            let match_mismatch = if seq1[i - 1] == seq2[j - 1] {
                self.match_score
            } else {
                self.mismatch_score
            };

            if current_score == score_matrix[i - 1][j - 1] + match_mismatch {
                align1.push(seq1[i - 1]);
                align2.push(seq2[j - 1]);
                i -= 1;
                j -= 1;
            } else if current_score == score_matrix[i - 1][j] + self.gap_penalty {
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
    for j in 0..=m {
        dp[j] = j;
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
}
