use exact_diagonalization::states::iterator::PeriodicPairIterator;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use exact_diagonalization::error::Error;

struct State {
    num: usize,
    length: usize,
}

struct StateIterator {
    num: usize,
    temp: usize,
    idx: usize,
    length: usize,
}

impl IntoIterator for State {
    type Item = (usize, usize);
    type IntoIter = StateIterator;

    fn into_iter(self) -> Self::IntoIter {
        StateIterator {
            num: (self.num / 2) + ((self.num % 2) << (self.length - 1)),
            temp: self.num % 2,
            idx: 0,
            length: self.length,
        }
    }
}

impl Iterator for StateIterator {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<(usize, usize)> {
        if self.idx == self.length {
            return None;
        } else {
            let result = (self.temp, self.num % 2);
            self.temp = self.num % 2;
            self.num = self.num / 2;
            self.idx += 1;
            return Some(result);
        }
    }
}

fn use_buildin(num: State) -> usize {
    // num = 1231412141232312
    // time:   [37.248 ns 38.727 ns 40.440 ns]
    // change: [+21.937% +29.110% +39.245%] (p = 0.00 < 0.05)
    let mut count = 0;
    for (_idx, (i, j)) in num.into_iter().enumerate() {
        count += i & j;
    }
    return count;
}

fn bench_buildin(c: &mut Criterion) {
    c.bench_function("use iterator", |b| {
        b.iter(|| {
            use_buildin(black_box(State {
                num: 1231412141232312,
                length: 51,
            }))
        })
    });
}

fn use_own_pkg(num: usize, length: usize) -> Result<usize, Error> {
    let mut count = 0;
    for (si, sj) in PeriodicPairIterator::new(num, length)? {
        count += si & sj;
    }
    return Ok(count);
}

fn bench_own_pkg(c: &mut Criterion) {
    c.bench_function("use iterator", |b| {
        b.iter(|| use_own_pkg(black_box(1231412141232312), black_box(51)))
    });
}

criterion_group!(benches, bench_buildin, bench_own_pkg);
criterion_main!(benches);
