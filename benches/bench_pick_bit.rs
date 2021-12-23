use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn use_and(num: usize, idx: usize) -> usize {
    // 8 digit : [1.1994 ns 1.2547 ns 1.3142 ns]
    // 20 digit : [1.1518 ns 1.2041 ns 1.2627 ns]
    return num & (1 << idx);
}

fn use_division(num: usize, idx: usize) -> usize {
    // 8 diti : [1.2266 ns 1.2803 ns 1.3407 ns]
    // 20 digit : [1.3753 ns 1.4695 ns 1.5802 ns]
    return (num / (1 << idx)) % 2;
}

fn bench_and(c: &mut Criterion) {
    c.bench_function("use and", |b| {
        b.iter(|| use_and(black_box(12412112541312), 20))
    });
}

fn bench_division(c: &mut Criterion) {
    c.bench_function("use division", |b| {
        b.iter(|| use_division(black_box(12412112541312), 20))
    });
}

criterion_group!(benches, bench_and, bench_division);
criterion_main!(benches);
