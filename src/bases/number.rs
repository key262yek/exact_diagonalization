use crate::{prelude::*, states::bit_fns::sum_bit};


pub type BasisN = BasisGenerator<EigenNumber>;
impl BasisN{
    pub fn new(v : EigenNumber, length : usize) -> Self{
        Self{
            length,
            value : Box::new(v),
        }
    }

    pub fn build(&self) -> Result<(Vec<NumberState>, FnvHashMap<RepNum, usize>), Error>{
        let num = self.value.total_number();
        let length = self.length;

        if num > length {
            return Err(Error::make_error_syntax(ErrorCode::OverFlow));
        }

        let max_state = 1 << length;
        let egn_v = &self.value;
        let mut basis : Vec<NumberState> = Vec::with_capacity(binomial(length, num));
        let mut indices : FnvHashMap<RepNum, usize> = FnvHashMap::default();
        let mut idx = 0;

        for n in 0..max_state {
            let state = SimpleState{rep : n, length};
            if state.bit_sum() == num {
                basis.push(NumberState::new(&state));
                indices.insert(Representation(**egn_v, n), idx);

                idx += 1;
            }
        }

        return Ok((basis, indices));
    }

    pub fn build_light(&self) -> Result<(Vec<(usize, usize)>, FnvHashMap<usize, usize>), Error>{
        let num = self.value.total_number();
        let length = self.length;

        if num > length {
            return Err(Error::make_error_syntax(ErrorCode::OverFlow));
        }

        let max_state = 1 << length;
        let mut basis : Vec<(usize, usize)> = Vec::with_capacity(binomial(length, num));
        let mut indices : FnvHashMap<usize, usize> = FnvHashMap::default();
        let mut idx = 0;

        for n in 0..max_state {
            if sum_bit(n) == num {
                basis.push((n, length));
                indices.insert(n, idx);
                idx += 1;
            }
        }

        return Ok((basis, indices));
    }
}


#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_basis_n(){
        let length = 4;
        let rep = 0;

        let gen = BasisN::new(EigenNumber::new(rep), length);
        let (base, indices) = gen.build().unwrap();
        assert_eq!(base, vec![NumberState::new(&SimpleState{rep : 0, length})]);
        for (&k, &v) in indices.iter(){
            assert!(k.get_rep() == 0 && v == 0);
        }

        let rep = 2;
        let gen = BasisN::new(EigenNumber::new(rep), length);
        let (base, indices) = gen.build().unwrap();
        assert_eq!(base,
            vec![NumberState::new(&SimpleState{rep : 3, length}),
                 NumberState::new(&SimpleState{rep : 5, length}),
                 NumberState::new(&SimpleState{rep : 6, length}),
                 NumberState::new(&SimpleState{rep : 9, length}),
                 NumberState::new(&SimpleState{rep : 10, length}),
                 NumberState::new(&SimpleState{rep : 12, length})]);
        for (&k, &v) in indices.iter(){
            assert!((k.get_rep() == 3 && v == 0) || (k.get_rep() == 5 && v == 1) || (k.get_rep() == 6 && v == 2)
                    || (k.get_rep() == 9 && v == 3) || (k.get_rep() == 10 && v == 4) || (k.get_rep() == 12 && v == 5));
        }
    }

    #[test]
    fn test_basis_light_n(){
        let length = 4;
        let rep = 0;

        let gen = BasisN::new(EigenNumber::new(rep), length);
        let (base, indices) = gen.build_light().unwrap();
        assert_eq!(base, vec![(0, 4)]);
        for (&k, &v) in indices.iter(){
            assert!(k == 0 && v == 0);
        }

        let rep = 2;
        let gen = BasisN::new(EigenNumber::new(rep), length);
        let (base, indices) = gen.build_light().unwrap();
        assert_eq!(base, vec![(3, 4), (5, 4), (6, 4), (9, 4), (10, 4), (12, 4)]);
        for (&k, &v) in indices.iter(){
            assert!((k == 3 && v == 0) || (k == 5 && v == 1) || (k == 6 && v == 2)
                    || (k == 9 && v == 3) || (k == 10 && v == 4) || (k == 12 && v == 5));
        }
    }
}
