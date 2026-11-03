//! Centrality algorithms for biological networks (Degree, Betweenness, Closeness).

use super::graph::BioGraph;
use std::collections::VecDeque;

/// Calculates normalized degree centrality for all nodes in the graph.
///
/// Degree centrality is defined as:
/// $C_D(v) = \frac{\text{deg}(v)}{N - 1}$ for $N > 1$.
pub fn degree_centrality<N, E: Clone>(graph: &BioGraph<N, E>) -> Vec<f64> {
    let n = graph.node_count();
    if n <= 1 {
        return vec![0.0; n];
    }
    let denom = (n - 1) as f64;
    (0..n).map(|i| graph.degree(i) as f64 / denom).collect()
}

/// Calculates betweenness centrality for all nodes using Ulrik Brandes' algorithm.
///
/// Computes shortest path counts and path dependencies in $O(V \cdot E)$ time.
pub fn betweenness_centrality<N, E: Clone>(
    graph: &BioGraph<N, E>,
    normalized: bool,
) -> Vec<f64> {
    let n = graph.node_count();
    let mut cb = vec![0.0; n];

    if n <= 2 {
        return cb;
    }

    for s in 0..n {
        let mut stack = Vec::with_capacity(n);
        let mut predecessors: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut sigma = vec![0.0; n];
        sigma[s] = 1.0;
        let mut dist = vec![-1i32; n];
        dist[s] = 0;

        let mut queue = VecDeque::new();
        queue.push_back(s);

        while let Some(v) = queue.pop_front() {
            stack.push(v);
            let d_v = dist[v];

            for &(w, _) in &graph.adjacency[v] {
                // w found for the first time?
                if dist[w] < 0 {
                    dist[w] = d_v + 1;
                    queue.push_back(w);
                }
                // Shortest path to w via v?
                if dist[w] == d_v + 1 {
                    sigma[w] += sigma[v];
                    predecessors[w].push(v);
                }
            }
        }

        let mut delta = vec![0.0; n];
        // S returns vertices in order of non-increasing distance from s
        while let Some(w) = stack.pop() {
            for &v in &predecessors[w] {
                if sigma[w] > 0.0 {
                    delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w]);
                }
            }
            if w != s {
                cb[w] += delta[w];
            }
        }
    }

    // For undirected graphs, every path was counted twice
    if !graph.is_directed {
        for val in &mut cb {
            *val /= 2.0;
        }
    }

    if normalized {
        let scale = if graph.is_directed {
            1.0 / ((n - 1) * (n - 2)) as f64
        } else {
            2.0 / ((n - 1) * (n - 2)) as f64
        };
        for val in &mut cb {
            *val *= scale;
        }
    }

    cb
}

/// Calculates closeness centrality for each node in the graph.
///
/// Uses Wasserman & Faust formulation for possibly disconnected graphs:
/// $C_C(u) = \frac{(R(u) - 1)^2}{(N - 1) \sum_{v \in R(u)} d(u, v)}$
/// where $R(u)$ is the number of reachable nodes from $u$.
pub fn closeness_centrality<N, E: Clone>(graph: &BioGraph<N, E>) -> Vec<f64> {
    let n = graph.node_count();
    let mut cc = vec![0.0; n];

    if n <= 1 {
        return cc;
    }

    for s in 0..n {
        let mut dist = vec![-1i32; n];
        dist[s] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(s);

        let mut reachable = 1usize;
        let mut sum_dist = 0.0;

        while let Some(u) = queue.pop_front() {
            let d_u = dist[u];
            for &(v, _) in &graph.adjacency[u] {
                if dist[v] < 0 {
                    dist[v] = d_u + 1;
                    reachable += 1;
                    sum_dist += dist[v] as f64;
                    queue.push_back(v);
                }
            }
        }

        if reachable > 1 && sum_dist > 0.0 {
            // Wasserman & Faust normalized closeness
            let reachable_f = (reachable - 1) as f64;
            let total_f = (n - 1) as f64;
            cc[s] = (reachable_f / total_f) * (reachable_f / sum_dist);
        }
    }

    cc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_star_graph_centrality() {
        // Star graph: Node 0 is center, connected to 1, 2, 3, 4
        let mut g = BioGraph::new(false);
        for i in 0..5 {
            g.add_node(format!("N{}", i));
        }
        for i in 1..5 {
            g.add_edge(0, i, 1.0).unwrap();
        }

        let deg = degree_centrality(&g);
        assert_eq!(deg[0], 1.0); // 4 / 4
        assert_eq!(deg[1], 0.25); // 1 / 4

        let bet = betweenness_centrality(&g, true);
        assert_eq!(bet[0], 1.0); // Center node mediates all pairs
        assert_eq!(bet[1], 0.0);

        let clo = closeness_centrality(&g);
        assert_eq!(clo[0], 1.0); // 1 hop to all 4 neighbors
        assert!(clo[1] < clo[0]);
    }
}
