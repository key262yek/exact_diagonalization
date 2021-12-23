// bitvec을 이용한 iterator와 현재의 division을 이용하는 방법 사이 속도 비교

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use bitvec::prelude::*;
use exact_diagonalization::states::{SimpleState, State};

fn simplestate(s : SimpleState) -> usize{
    // Collecting 100 samples in estimated 5.0003 s (47M iterations)
    // time:   [105.35 ns 105.62 ns 105.88 ns]

    let mut sum = 0;
    for i in s.bit_iter(){
        sum += i;
    }
    return sum;
}

fn using_bitvec(n : usize) -> usize{
    // Collecting 100 samples in estimated 5.0009 s (13M iterations)
    // time:   [388.87 ns 389.22 ns 389.58 ns]
    let mut s : BitVec<Lsb0, usize> = BitVec::new();
    s.resize(64, false);
    s[..64].store(n);
    let mut sum = 0;
    for i in s{
        sum += i as usize;
    }

    return sum;
}

fn bench_bitvec(c: &mut Criterion) {
    c.bench_function("bitvec iterator", |b| {
        b.iter(|| {
            using_bitvec(black_box(1231412141232312))
        })
    });
}

fn bench_simplestate(c: &mut Criterion) {
    c.bench_function("simple state iterator", |b| {
        b.iter(|| simplestate(SimpleState{rep : 1231412141232312, length : 51}))
    });
}

criterion_group!(benches, bench_bitvec, bench_simplestate);
criterion_main!(benches);
