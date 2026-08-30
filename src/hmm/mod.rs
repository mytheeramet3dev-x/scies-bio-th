/// State kinds in a Profile HMM.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum HmmState {
    Match(usize),
    Insert(usize),
    Delete(usize),
    Start,
    End,
}

fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// A simple linear Profile HMM with emission and transition probabilities.
///
/// Architecture: Start -> M1/I0 -> M2/I1 -> ... -> Mk/Ik -> End
/// For each column i there are states: Match(i), Insert(i), Delete(i).
pub struct ProfileHmm {
    /// Number of match states (= length of the alignment profile).
    pub length: usize,
    /// Emission probabilities for Match states. match_emit[i][ACGT].
    pub match_emit: Vec<[f64; 4]>,
    /// Emission probabilities for Insert states (uniform 0.25 default).
    pub insert_emit: Vec<[f64; 4]>,
    /// Transition probabilities from each column's Match state.
    /// trans[i] = [p(M->M), p(M->I), p(M->D)]
    pub trans_match: Vec<[f64; 3]>,
    /// Transition probabilities from each column's Insert state.
    /// trans_insert[i] = [p(I->M), p(I->I)]
    pub trans_insert: Vec<[f64; 2]>,
}

impl ProfileHmm {
    /// Creates a new Profile HMM with given length and uniform priors.
    pub fn new(length: usize) -> Self {
        Self {
            length,
            match_emit: vec![[0.25; 4]; length],
            insert_emit: vec![[0.25; 4]; length + 1],
            trans_match: vec![[0.7, 0.1, 0.2]; length],
            trans_insert: vec![[0.8, 0.2]; length + 1],
        }
    }

    /// Runs the Viterbi algorithm on a DNA byte slice.
    /// Returns the (score, path of states) with the highest probability.
    pub fn viterbi(&self, seq: &[u8]) -> (f64, Vec<HmmState>) {
        let n = seq.len();
        let k = self.length;

        // dp[i][j] = log-prob of best path ending at seq[i] being emitted by Match(j)
        // We'll track only Match states for simplicity in this v0.2 implementation.
        let neg_inf = f64::NEG_INFINITY;
        let mut dp = vec![vec![neg_inf; k + 1]; n + 1];
        let mut back: Vec<Vec<Option<usize>>> = vec![vec![None; k + 1]; n + 1];

        dp[0][0] = 0.0; // Start state

        for i in 1..=n {
            let bi = base_index(seq[i - 1]);
            for j in 1..=k {
                // Transition from Match(j-1) -> Match(j)
                let emit_score = bi
                    .map(|idx| self.match_emit[j - 1][idx].ln())
                    .unwrap_or(neg_inf);

                let p_trans = self.trans_match[j - 1][0].ln(); // M -> M
                let candidate = dp[i - 1][j - 1] + p_trans + emit_score;
                if candidate > dp[i][j] {
                    dp[i][j] = candidate;
                    back[i][j] = Some(j - 1);
                }
            }
        }

        // Find best ending state
        let (best_j, best_score) = (1..=k)
            .map(|j| (j, dp[n][j]))
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap_or((0, neg_inf));

        // Traceback
        let mut path = Vec::new();
        let mut cur_i = n;
        let mut cur_j = best_j;
        while cur_i > 0 && cur_j > 0 {
            path.push(HmmState::Match(cur_j));
            if let Some(prev_j) = back[cur_i][cur_j] {
                cur_j = prev_j;
                cur_i -= 1;
            } else {
                break;
            }
        }
        path.push(HmmState::Start);
        path.reverse();

        (best_score, path)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hmm_viterbi_basic() {
        let mut hmm = ProfileHmm::new(4);
        // Tune match emissions to strongly prefer ACGT
        hmm.match_emit[0] = [0.97, 0.01, 0.01, 0.01]; // A
        hmm.match_emit[1] = [0.01, 0.97, 0.01, 0.01]; // C
        hmm.match_emit[2] = [0.01, 0.01, 0.97, 0.01]; // G
        hmm.match_emit[3] = [0.01, 0.01, 0.01, 0.97]; // T

        let (score, path) = hmm.viterbi(b"ACGT");
        assert!(score > f64::NEG_INFINITY, "Should find a valid path");
        assert!(path.len() > 0);
        assert_eq!(path[0], HmmState::Start);
    }

    #[test]
    fn test_hmm_viterbi_mismatch_lower_score() {
        let mut hmm = ProfileHmm::new(4);
        hmm.match_emit[0] = [0.97, 0.01, 0.01, 0.01]; // A
        hmm.match_emit[1] = [0.01, 0.97, 0.01, 0.01]; // C
        hmm.match_emit[2] = [0.01, 0.01, 0.97, 0.01]; // G
        hmm.match_emit[3] = [0.01, 0.01, 0.01, 0.97]; // T

        let (score_match, _) = hmm.viterbi(b"ACGT");
        let (score_mismatch, _) = hmm.viterbi(b"TTTT");
        assert!(score_match > score_mismatch, "Perfect match should score higher");
    }
}
