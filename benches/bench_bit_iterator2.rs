use criterion::{black_box, criterion_group, criterion_main, Criterion};

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

fn use_iterator(num: State) -> usize {
    // num = 1231412141232312
    // time:   [32.555 ns 32.798 ns 33.059 ns]
    // change: [-60.004% -59.100% -58.264%] (p = 0.00 < 0.05)
    let mut count = 0;
    for (i, j) in num.into_iter() {
        count += i & j;
    }
    return count;
}

fn use_pick_bit(num: usize, length: usize) -> usize {
    // num = 1231412141232312, length = 51
    // time:   [58.217 ns 58.629 ns 59.217 ns]
    // change: [-4.2014% -2.3838% -0.5275%] (p = 0.01 < 0.05)
    fn pick_bit(num: usize, idx: usize) -> usize {
        return num & (1 << idx);
    }

    let mut count = 0;
    let mut temp = pick_bit(num, 0);
    let mut temp2: usize;
    for i in 1..length {
        temp2 = pick_bit(num, i);
        count += temp & temp2;
        temp = temp2;
    }
    count += pick_bit(num, 0) & pick_bit(num, length - 1);

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

criterion_group!(benches, bench_iterator, bench_pick_bit);
criterion_main!(benches);
