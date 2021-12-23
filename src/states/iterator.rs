use crate::prelude::*;

#[derive(Copy, Clone, Debug, Eq, PartialEq, PartialOrd)]
pub struct BitIterator {
    num: usize,
    idx: usize,
    length: usize,
}

impl BitIterator {
    #[allow(dead_code)]
    pub fn new(num: usize, length: usize) -> Result<Self, Error> {
        if num >= (1 << length) {
            return Err(Error::make_error_syntax(ErrorCode::OverFlow));
        }
        Ok(Self {
            num,
            idx: 0,
            length,
        })
    }
}

impl Iterator for BitIterator {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx == self.length {
            return None;
        } else {
            let temp = self.num % 2;
            self.idx += 1;
            self.num = self.num >> 1;
            return Some(temp);
        }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, PartialOrd)]
pub struct PairIterator {
    num: usize,
    temp: usize,
    idx: usize,
    length: usize,
}

impl PairIterator {
    #[allow(dead_code)]
    pub fn new(num: usize, length: usize) -> Result<Self, Error> {
        if num >= (1 << length) {
            return Err(Error::make_error_syntax(ErrorCode::OverFlow));
        }
        Ok(Self {
            num: num >> 1,
            temp: num % 2,
            idx: 0,
            length: length - 1,
        })
    }
}

impl Iterator for PairIterator {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx == self.length {
            return None;
        } else {
            let result = (self.temp, self.num % 2);
            self.temp = self.num % 2;
            self.num = self.num >> 1;
            self.idx += 1;
            return Some(result);
        }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, PartialOrd)]
pub struct PeriodicPairIterator {
    num: usize,
    temp: usize,
    idx: usize,
    length: usize,
}

impl PeriodicPairIterator {
    #[allow(dead_code)]
    pub fn new(num: usize, length: usize) -> Result<Self, Error> {
        if num >= (1 << length) {
            return Err(Error::make_error_syntax(ErrorCode::OverFlow));
        }
        Ok(Self {
            num: (num >> 1) + ((num % 2) << (length - 1)),
            temp: num % 2,
            idx: 0,
            length,
        })
    }
}

impl Iterator for PeriodicPairIterator {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx == self.length {
            return None;
        } else {
            let result = (self.temp, self.num % 2);
            self.temp = self.num % 2;
            self.num = self.num >> 1;
            self.idx += 1;
            return Some(result);
        }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, PartialOrd)]
pub struct PeriodicPairEnumerator {
    num: usize,
    temp: usize,
    idx: usize,
    length: usize,
}

impl PeriodicPairEnumerator {
    #[allow(dead_code)]
    pub fn new(num: usize, length: usize) -> Result<Self, Error> {
        if num >= (1 << length) {
            return Err(Error::make_error_syntax(ErrorCode::OverFlow));
        }
        Ok(Self {
            num: (num >> 1) + ((num % 2) << (length - 1)),
            temp: num % 2,
            idx: 0,
            length,
        })
    }
}

impl Iterator for PeriodicPairEnumerator {
    type Item = ((usize, usize), (usize, usize));

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx == self.length {
            return None;
        } else if self.idx == self.length - 1{
            let result = ((self.idx, self.temp), (0, self.num % 2));
            self.idx += 1;
            return Some(result);
        } else {
            let result = ((self.idx, self.temp), (self.idx + 1, self.num % 2));
            self.temp = self.num % 2;
            self.num = self.num >> 1;
            self.idx += 1;
            return Some(result);
        }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, PartialOrd)]
pub struct PeriodicDistancedPairIterator{
    num1 : usize,
    num2 : usize,
    idx : usize,
    length : usize,
    dist : usize
}

impl PeriodicDistancedPairIterator{
    pub fn new(num : usize, length : usize, dist : usize) -> Result<Self, Error>{
        if num >= (1 << length) || dist >= length {
            return Err(Error::make_error_syntax(ErrorCode::OverFlow));
        } else {
            Ok(Self{
                num1 : num,
                num2 : ((num >> dist) + ((num % (1 << dist)) << (length - dist))),
                idx : 0,
                length,
                dist,
            })
        }
    }
}

impl Iterator for PeriodicDistancedPairIterator{
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item>{
        if self.idx == self.length{
            return None;
        } else {
            let result = (self.num1 % 2, self.num2 % 2);
            self.num1 = self.num1 >> 1;
            self.num2 = self.num2 >> 1;
            self.idx += 1;
            return Some(result);
        }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, PartialOrd)]
pub struct PeriodicDistancedPairEnumerator{
    num1 : usize,
    num2 : usize,
    idx : usize,
    length : usize,
    dist : usize
}

impl PeriodicDistancedPairEnumerator{
    pub fn new(num : usize, length : usize, dist : usize) -> Result<Self, Error>{
        if num >= (1 << length) || dist >= length {
            return Err(Error::make_error_syntax(ErrorCode::OverFlow));
        } else {
            Ok(Self{
                num1 : num,
                num2 : ((num >> dist) + ((num % (1 << dist)) << (length - dist))),
                idx : 0,
                length,
                dist,
            })
        }
    }
}

impl Iterator for PeriodicDistancedPairEnumerator{
    type Item = ((usize, usize), (usize, usize));

    fn next(&mut self) -> Option<Self::Item>{
        if self.idx == self.length{
            return None;
        } else {
            let result = ((self.idx, self.num1 % 2), ((self.idx + self.dist) % self.length, self.num2 % 2));
            self.num1 = self.num1 >> 1;
            self.num2 = self.num2 >> 1;
            self.idx += 1;
            return Some(result);
        }
    }
}


#[derive(Copy, Clone, Debug, Eq, PartialEq, PartialOrd)]
pub struct CycleIterator{
    start : usize,
    num : usize,
    idx : usize,
    length : usize
}

impl CycleIterator{
    #[allow(dead_code)]
    pub fn new(num : usize, length : usize) -> Result<Self, Error>{
        if num >= (1 << length) {
            return Err(Error::make_error_syntax(ErrorCode::OverFlow));
        }
        Ok(Self{
            start : num,
            num,
            idx : 0,
            length,
        })
    }
}

impl Iterator for CycleIterator{
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item>{
        if self.idx == 0 {
            self.idx += 1;
            self.num = cyclic_move_unsafe(self.num, self.length);
            return Some(self.start);
        }
        else if self.idx == self.length || self.start == self.num{
            return None;
        } else {
            let result = self.num;
            self.num = cyclic_move_unsafe(self.num, self.length);
            self.idx += 1;
            return Some(result);
        }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, PartialOrd)]
pub struct CommenIterator{
    current : usize,
    period : usize,
    length : usize,
}

impl CommenIterator{
    pub fn new(period : usize, length : usize) -> Self{
        Self{
            current : 0,
            period,
            length,
        }
    }
}

impl Iterator for CommenIterator{
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item>{
        let mut num = self.current;
        while num < self.length {
            if (num * self.period) % self.length == 0{
                self.current = num + 1;
                return Some(num);
            }
            num += 1;
        }
        return None;
    }
}

#[cfg(test)]
mod test {
    use crate::states::{SimpleState, State};
    use super::*;

    #[test]
    fn test_simple_iter(){
        let state = SimpleState::new(10, 5);

        let mut x = state.bit_iter();
        assert_eq!(x.next(), Some(0));
        assert_eq!(x.next(), Some(1));
        assert_eq!(x.next(), Some(0));
        assert_eq!(x.next(), Some(1));
        assert_eq!(x.next(), Some(0));
        assert_eq!(x.next(), None);

        let mut x = state.pair_iter();
        assert_eq!(x.next(), Some((0, 1)));
        assert_eq!(x.next(), Some((1, 0)));
        assert_eq!(x.next(), Some((0, 1)));
        assert_eq!(x.next(), Some((1, 0)));
        assert_eq!(x.next(), None);

        let mut x = state.periodic_pair_iter();
        assert_eq!(x.next(), Some((0, 1)));
        assert_eq!(x.next(), Some((1, 0)));
        assert_eq!(x.next(), Some((0, 1)));
        assert_eq!(x.next(), Some((1, 0)));
        assert_eq!(x.next(), Some((0, 0)));
        assert_eq!(x.next(), None);

        let mut x = state.periodic_pair_enumerate();
        assert_eq!(x.next(), Some(((0, 0), (1, 1))));
        assert_eq!(x.next(), Some(((1, 1), (2, 0))));
        assert_eq!(x.next(), Some(((2, 0), (3, 1))));
        assert_eq!(x.next(), Some(((3, 1), (4, 0))));
        assert_eq!(x.next(), Some(((4, 0), (0, 0))));
        assert_eq!(x.next(), None);

        let mut x = state.cycle_iter();
        assert_eq!(x.next(), Some(10));
        assert_eq!(x.next(), Some(5));
        assert_eq!(x.next(), Some(18));
        assert_eq!(x.next(), Some(9));
        assert_eq!(x.next(), Some(20));
        assert_eq!(x.next(), None);
    }

    #[test]
    fn test_bit_iter() -> Result<(), Error> {
        assert_eq!(
            BitIterator::new(10, 3),
            Err(Error::make_error_syntax(ErrorCode::OverFlow))
        );

        let mut bit_iterator = BitIterator::new(10, 4)?;
        assert_eq!(bit_iterator.next(), Some(0));
        assert_eq!(bit_iterator.next(), Some(1));
        assert_eq!(bit_iterator.next(), Some(0));
        assert_eq!(bit_iterator.next(), Some(1));
        assert_eq!(bit_iterator.next(), None);

        let mut bit_iterator = BitIterator::new(10, 5)?;
        assert_eq!(bit_iterator.next(), Some(0));
        assert_eq!(bit_iterator.next(), Some(1));
        assert_eq!(bit_iterator.next(), Some(0));
        assert_eq!(bit_iterator.next(), Some(1));
        assert_eq!(bit_iterator.next(), Some(0));
        assert_eq!(bit_iterator.next(), None);

        return Ok(());
    }

    #[test]
    fn test_pair_iter() -> Result<(), Error> {
        assert_eq!(
            PairIterator::new(10, 3),
            Err(Error::make_error_syntax(ErrorCode::OverFlow))
        );

        let mut pair_iterator = PairIterator::new(10, 4)?;
        assert_eq!(pair_iterator.next(), Some((0, 1)));
        assert_eq!(pair_iterator.next(), Some((1, 0)));
        assert_eq!(pair_iterator.next(), Some((0, 1)));
        assert_eq!(pair_iterator.next(), None);

        let mut pair_iterator = PairIterator::new(10, 5)?;
        assert_eq!(pair_iterator.next(), Some((0, 1)));
        assert_eq!(pair_iterator.next(), Some((1, 0)));
        assert_eq!(pair_iterator.next(), Some((0, 1)));
        assert_eq!(pair_iterator.next(), Some((1, 0)));
        assert_eq!(pair_iterator.next(), None);

        return Ok(());
    }

    #[test]
    fn test_periodic_pair_iter() -> Result<(), Error> {
        assert_eq!(
            PeriodicPairIterator::new(10, 3),
            Err(Error::make_error_syntax(ErrorCode::OverFlow))
        );

        let mut period_iter = PeriodicPairIterator::new(10, 4)?;
        assert_eq!(period_iter.next(), Some((0, 1)));
        assert_eq!(period_iter.next(), Some((1, 0)));
        assert_eq!(period_iter.next(), Some((0, 1)));
        assert_eq!(period_iter.next(), Some((1, 0)));
        assert_eq!(period_iter.next(), None);

        let mut period_iter = PeriodicPairIterator::new(10, 5)?;
        assert_eq!(period_iter.next(), Some((0, 1)));
        assert_eq!(period_iter.next(), Some((1, 0)));
        assert_eq!(period_iter.next(), Some((0, 1)));
        assert_eq!(period_iter.next(), Some((1, 0)));
        assert_eq!(period_iter.next(), Some((0, 0)));
        assert_eq!(period_iter.next(), None);

        return Ok(());
    }

    #[test]
    fn test_periodic_pair_enumerate() -> Result<(), Error> {
        assert_eq!(
            PeriodicPairEnumerator::new(10, 3),
            Err(Error::make_error_syntax(ErrorCode::OverFlow))
        );

        let mut period_iter = PeriodicPairEnumerator::new(10, 4)?;
        assert_eq!(period_iter.next(), Some(((0, 0), (1, 1))));
        assert_eq!(period_iter.next(), Some(((1, 1), (2, 0))));
        assert_eq!(period_iter.next(), Some(((2, 0), (3, 1))));
        assert_eq!(period_iter.next(), Some(((3, 1), (0, 0))));
        assert_eq!(period_iter.next(), None);

        let mut period_iter = PeriodicPairEnumerator::new(10, 5)?;
        assert_eq!(period_iter.next(), Some(((0, 0), (1, 1))));
        assert_eq!(period_iter.next(), Some(((1, 1), (2, 0))));
        assert_eq!(period_iter.next(), Some(((2, 0), (3, 1))));
        assert_eq!(period_iter.next(), Some(((3, 1), (4, 0))));
        assert_eq!(period_iter.next(), Some(((4, 0), (0, 0))));
        assert_eq!(period_iter.next(), None);

        return Ok(());
    }

    #[test]
    fn test_cycle_iter() -> Result<(), Error> {
        assert_eq!(
            CycleIterator::new(10, 3),
            Err(Error::make_error_syntax(ErrorCode::OverFlow))
        );

        let mut cycle_it = CycleIterator::new(10, 4)?;
        assert_eq!(cycle_it.next(), Some(10));
        assert_eq!(cycle_it.next(), Some(5));
        assert_eq!(cycle_it.next(), None);

        let mut cycle_it = CycleIterator::new(10, 5)?;
        assert_eq!(cycle_it.next(), Some(10));
        assert_eq!(cycle_it.next(), Some(5));
        assert_eq!(cycle_it.next(), Some(18));
        assert_eq!(cycle_it.next(), Some(9));
        assert_eq!(cycle_it.next(), Some(20));
        assert_eq!(cycle_it.next(), None);

        return Ok(());
    }

    #[test]
    fn test_commensurability_iterator(){
        let mut s = CommenIterator::new(6, 12);
        assert_eq!(s.next(), Some(0));
        assert_eq!(s.next(), Some(2));
        assert_eq!(s.next(), Some(4));
        assert_eq!(s.next(), Some(6));
        assert_eq!(s.next(), Some(8));
        assert_eq!(s.next(), Some(10));
        assert_eq!(s.next(), None);

        let mut s = CommenIterator::new(4, 12);
        assert_eq!(s.next(), Some(0));
        assert_eq!(s.next(), Some(3));
        assert_eq!(s.next(), Some(6));
        assert_eq!(s.next(), Some(9));
        assert_eq!(s.next(), None);

    }
}
