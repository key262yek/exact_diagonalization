use criterion::{black_box, criterion_group, criterion_main, Criterion};
use exact_diagonalization::states::{SimpleState, bit_fns::bit_flip};
use exact_diagonalization::error::Error;

fn normal_number(num : usize, length : usize) -> Result<(), Error>{
    // Collecting 100 samples in estimated 5.0002 s (41M iterations)
    // time:   [120.21 ns 122.93 ns 126.38 ns]
    // 6 (6.00%) high mild
    // 7 (7.00%) high severe
    let mut temp  = num;
    for i in 0..length-1{
        temp = bit_flip(temp, length, i, i + 1)?;
    }

    Ok(())
}

fn simplestate(s : SimpleState) -> Result<(), Error>{
    // Collecting 100 samples in estimated 5.0015 s (14M iterations)
    // time:   [352.32 ns 360.95 ns 371.42 ns]
    // 5 (5.00%) high mild
    // 8 (8.00%) high severe
    let mut temp = s;
    for i in 0..s.length - 1{
        temp = temp.bit_flip(i, i + 1)?;
    }
    Ok(())
}

fn simplestate2(s : SimpleState) -> Result<(), Error>{
    // Collecting 100 samples in estimated 5.0003 s (43M iterations)
    // time:   [114.41 ns 117.08 ns 120.50 ns]
    // 3 (3.00%) high mild
    // 6 (6.00%) high severe
    let mut temp = s;
    for i in 0..s.length - 1{
        temp.rep = bit_flip(temp.rep, temp.length, i, i + 1)?;
    }
    Ok(())
}

fn simplestate3(s : SimpleState){
    // Collecting 100 samples in estimated 5.0003 s (46M iterations)
    // time:   [107.76 ns 111.98 ns 117.33 ns]
    // 3 (3.00%) high mild
    // 11 (11.00%) high severe
    let mut temp = s;
    for i in 0..s.length - 1{
        temp.bit_flip_mut(i, i + 1);
    }
}

fn bench_normal(c: &mut Criterion) {
    c.bench_function("normal number allocation", |b| {
        b.iter(|| {
            normal_number(black_box(1231412141232312), black_box(51))
        })
    });
}

fn bench_simplestate(c: &mut Criterion) {
    c.bench_function("simple state allocation", |b| {
        b.iter(|| simplestate(SimpleState{rep : 1231412141232312, length : 51}))
    });
}

fn bench_simplestate2(c: &mut Criterion) {
    c.bench_function("simple state allocation", |b| {
        b.iter(|| simplestate2(SimpleState{rep : 1231412141232312, length : 51}))
    });
}

fn bench_simplestate3(c: &mut Criterion) {
    c.bench_function("simple state allocation", |b| {
        b.iter(|| simplestate3(SimpleState{rep : 1231412141232312, length : 51}))
    });
}

criterion_group!(benches, bench_normal, bench_simplestate, bench_simplestate2, bench_simplestate3);
criterion_main!(benches);
