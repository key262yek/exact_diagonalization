use criterion::{black_box, criterion_group, criterion_main, Criterion};
use exact_diagonalization::states::SimpleState;

fn normal_number(num : usize) -> usize{
    // Collecting 100 samples in estimated 5.0000 s (276M iterations)
    // time:   [17.673 ns 18.174 ns 18.770 ns]
    // change: [-86.793% -86.392% -85.918%] (p = 0.00 < 0.05)
    let mut temp  = num;
    let mut sum = 0;
    while temp != 0{
        sum += temp % 2;
        temp /= 2;
    }

    sum
}

fn simplestate(s : SimpleState) -> usize{
    // Collecting 100 samples in estimated 5.0000 s (279M iterations)
    // time:   [17.921 ns 18.563 ns 19.377 ns]
    // change: [-94.790% -94.573% -94.349%] (p = 0.00 < 0.05)
    s.total_number()
}

fn bench_normal(c: &mut Criterion) {
    c.bench_function("normal number allocation", |b| {
        b.iter(|| {
            normal_number(black_box(1231412141232312))
        })
    });
}

fn bench_simplestate(c: &mut Criterion) {
    c.bench_function("simple state allocation", |b| {
        b.iter(|| simplestate(SimpleState{rep : 1231412141232312, length : 51}))
    });
}

criterion_group!(benches, bench_normal, bench_simplestate);
criterion_main!(benches);
