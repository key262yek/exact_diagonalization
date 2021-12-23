use crate::prelude::*;

impl SimpleState{
    pub fn bit_flip(&self, i : usize, j : usize) -> Result<Self, Error>{
        bit_flip(self.rep, self.length, i, j).map(|x| Self{rep : x, length : self.length})
    }

    pub fn bit_flip_mut(&mut self, i : usize, j : usize){
        match bit_flip(self.rep, self.length, i, j){
            Ok(x) => {self.rep = x;
            },
            Err(e) => {panic!("{}", e);
            },
        }
    }

    pub fn pick_bit(&self, idx : usize) -> Result<usize, Error>{
        pick_bit(self.rep, self.length, idx)
    }

    pub fn cyclic_move(&self) -> Self{
        let rep = (self.rep >> 1) + ((self.rep % 2) << (self.length - 1));
        Self{
            rep,
            length : self.length
        }
    }

    pub fn cyclic_move_mut(&mut self){
        let n = self.rep;
        self.rep = (n >> 1) + ((n % 2) << (self.length - 1));
    }
}

pub fn bit_flip(num: usize, length: usize, i: usize, j: usize) -> Result<usize, Error> {
    if length <= i || length <= j || i == j {
        return Err(Error::make_error_syntax(ErrorCode::InvalidBitIndex));
    }

    return Ok(num ^ ((1 << i) + (1 << j)));
}

pub fn bit_flip_unsafe(num: usize, i: usize, j: usize) -> usize{
    return num ^ ((1 << i) + (1 << j));
}

pub fn pick_bit(num: usize, length: usize, idx: usize) -> Result<usize, Error> {
    if length <= idx {
        return Err(Error::make_error_syntax(ErrorCode::InvalidBitIndex));
    }

    return Ok(((num >> idx) % 2) as usize);
}

pub fn sum_bit(num: usize) -> usize {
    let mut temp = num;
    let mut count = 0;
    while temp > 0 {
        count += temp % 2;
        temp = temp >> 1;
    }

    return count;
}

pub fn cyclic_move(num: usize, length: usize) -> Result<usize, Error> {
    // cyclic move like 10010 => 01001 => 10100 => 01010 => 00101 => 10010
    if num >= (1 << length) {
        return Err(Error::make_error_syntax(ErrorCode::OverFlow));
    }

    return Ok((num >> 1) + ((num % 2) << (length - 1)));
}

pub fn cyclic_move_unsafe(num : usize, length : usize) -> usize{
    return (num >> 1) + ((num % 2) << (length - 1));
}

pub fn period(num: usize, length: usize) -> Result<usize, Error> {
    if num >= (1 << length) {
        return Err(Error::make_error_syntax(ErrorCode::OverFlow));
    }

    let mut temp = num;
    let mut count = 0;

    loop {
        temp = cyclic_move(temp, length)?;
        count += 1;
        if temp == num {
            return Ok(count);
        }
    }
}

pub fn period_unsafe(num : usize, length : usize) -> usize{
    let mut temp = num;
    let mut count = 0;

    loop {
        temp = cyclic_move_unsafe(temp, length);
        count += 1;
        if temp == num {
            return count;
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_bit_fn() {
        assert_eq!(pick_bit(10, 4, 0), Ok(0));
        assert_eq!(pick_bit(10, 4, 1), Ok(1));
        assert_eq!(pick_bit(10, 4, 2), Ok(0));
        assert_eq!(pick_bit(10, 4, 3), Ok(1));
        assert_eq!(
            pick_bit(10, 4, 4),
            Err(Error::make_error_syntax(ErrorCode::OverFlow))
        );
        assert_eq!(pick_bit(10, 5, 4), Ok(0));

        assert_eq!(sum_bit(4), 1);
        assert_eq!(sum_bit(10), 2);
        assert_eq!(sum_bit(7), 3);
        assert_eq!(sum_bit(31), 5);

        assert_eq!(cyclic_move(4, 3), Ok(2));
        assert_eq!(cyclic_move(2, 3), Ok(1));
        assert_eq!(cyclic_move(1, 3), Ok(4));

        assert_eq!(cyclic_move(18, 5), Ok(9));
        assert_eq!(cyclic_move(9, 5), Ok(20));
        assert_eq!(cyclic_move(20, 5), Ok(10));
        assert_eq!(cyclic_move(10, 5), Ok(5));
        assert_eq!(cyclic_move(5, 5), Ok(18));
        assert_eq!(
            cyclic_move(18, 3),
            Err(Error::make_error_syntax(ErrorCode::OverFlow))
        );

        assert_eq!(period(4, 3), Ok(3));
        assert_eq!(period(10, 4), Ok(2));
        assert_eq!(period(10, 5), Ok(5));
        assert_eq!(period(36, 6), Ok(3));
        assert_eq!(period(7, 3), Ok(1));
        assert_eq!(
            period(18, 3),
            Err(Error::make_error_syntax(ErrorCode::OverFlow))
        );
    }
}
