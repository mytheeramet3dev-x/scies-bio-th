use criterion::{black_box, criterion_group, criterion_main, Criterion};
use scies_bio_th::align::{levenshtein_distance, NeedlemanWunsch, SmithWaterman};

pub fn align_benchmark(c: &mut Criterion) {
    let seq1 = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
    let seq2 = b"ACGTACGTTCGTACGTACGTCCGTACGTACGA";

    let nw = NeedlemanWunsch::new(1, -1, -2);
    c.bench_function("needleman_wunsch_32bp", |b| {
        b.iter(|| nw.align(black_box(seq1), black_box(seq2)))
    });

    let sw = SmithWaterman::new(2, -1, -2);
    c.bench_function("smith_waterman_32bp", |b| {
        b.iter(|| sw.align(black_box(seq1), black_box(seq2)))
    });

    c.bench_function("levenshtein_32bp", |b| {
        b.iter(|| levenshtein_distance(black_box(seq1), black_box(seq2)))
    });
}

criterion_group!(benches, align_benchmark);
criterion_main!(benches);
