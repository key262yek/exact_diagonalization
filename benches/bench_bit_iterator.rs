use criterion::{black_box, criterion_group, criterion_main, Criterion};

struct State {
    num: usize,
    length: usize,
}

struct StateIterator {
    num: usize,
    idx: usize,
    length: usize,
}

impl IntoIterator for State {
    type Item = usize;
    type IntoIter = StateIterator;

    fn into_iter(self) -> Self::IntoIter {
        StateIterator {
            num: self.num,
            idx: 0,
            length: self.length,
        }
    }
}

impl Iterator for StateIterator {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        if self.idx == self.length {
            return None;
        } else {
            let temp = self.num % 2;
            self.idx += 1;
            self.num = self.num / 2;
            return Some(temp);
        }
    }
}

fn use_iterator(num: State) -> usize {
    // num = 1231412141232312
    // Collecting 100 samples in estimated 5.0000 s (321M iterations)
    // time:   [15.614 ns 15.651 ns 15.689 ns]
    // change: [-59.366% -57.797% -56.259%] (p = 0.00 < 0.05)
    let mut count = 0;
    for i in num.into_iter() {
        count += i;
    }
    return count;
}

fn use_pick_bit(num: usize, length: usize) -> usize {
    // num = 1231412141232312, length = 51
    // Collecting 100 samples in estimated 5.0001 s (190M iterations)
    // time:   [26.124 ns 26.161 ns 26.196 ns]
    // change: [-57.979% -57.375% -56.824%] (p = 0.00 < 0.05)
    fn pick_bit(num: usize, idx: usize) -> usize {
        return num & (1 << idx);
    }

    let mut count = 0;
    for i in 0..length {
        count += pick_bit(num, i);
    }

    return count;
}

fn use_bit_and(num: usize, length: usize) -> usize {
    // num = 1231412141232312, length = 51
    // Collecting 100 samples in estimated 5.0000 s (366M iterations)
    // time:   [13.543 ns 13.563 ns 13.580 ns]

    let mut count = 0;
    let mut mask = 1;
    for _ in 0..length {
        count += num & mask;
        mask = mask << 1;
    }

    return count;
}

fn bench_iterator(c: &mut Criterion) {
    c.bench_function("use iterator", |b| {
        b.iter(|| {
            use_iterator(black_box(State {
                num: 1231412141232312,
                length: 51,
            }))
        })
    });
}

fn bench_pick_bit(c: &mut Criterion) {
    c.bench_function("use pick_bit", |b| {
        b.iter(|| use_pick_bit(black_box(1231412141232312), black_box(51)))
    });
}

fn bench_bit_and(c: &mut Criterion) {
    c.bench_function("use pick_bit", |b| {
        b.iter(|| use_bit_and(black_box(1231412141232312), black_box(51)))
    });
}

criterion_group!(benches, bench_iterator, bench_pick_bit, bench_bit_and);
criterion_main!(benches);
