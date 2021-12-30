use crate::{prelude::*, states::bit_fns::{sum_bit, is_rep}};

pub type BasisNK = BasisGenerator<EigenNumMomentum>;

impl BasisNK{
    pub fn new(v : EigenNumMomentum, length : usize) -> Self{
        Self{
            length,
            value : Box::new(v),
        }
    }

    pub fn check_commensurability(&self, period : usize) -> bool{
        self.value.check_commensurability(period, self.length)
    }

    pub fn build(&self) -> Result<(Vec<NumMomentumState>, FnvHashMap<usize, (usize, usize)>), Error>{
        let max_state = 1 << self.length;
        let mut basis : Vec<NumMomentumState> = Vec::new();
        let mut indices : FnvHashMap<usize, (usize, usize)> = FnvHashMap::default();
        let mut idx = 0;

        for n in 0..max_state{
            let state = SimpleState::new(n, self.length);
            let m = state.total_number();

            if m != self.value.total_number() {
                continue;
            }

            if let Some(nkstate) = NumMomentumState::new(state, self.value.wave_number()){
                for (i, state) in nkstate.state.states.iter().enumerate(){
                    let num = state.rep;
                    indices.insert(num, (idx, i));
                }

                basis.push(nkstate);
                idx += 1;
            }
        }

        if basis.len() == 0{
            return Err(Error::make_error_syntax(ErrorCode::InvalidConfiguration));
        }
        return Ok((basis, indices));
    }

    pub fn build_light(&self) -> Result<(Vec<(usize, usize)>, FnvHashMap<usize, (usize, usize)>), Error>{
        let length = self.length;
        let max_state = 1 << length;
        let mut basis : Vec<(usize, usize)> = Vec::new();
        let mut indices : FnvHashMap<usize, (usize, usize)> = FnvHashMap::default();
        let mut idx = 0;

        for n in 0..max_state{
            let m = sum_bit(n);
            let p = period_unsafe(n, length);

            if m != self.value.total_number()
                || !self.check_commensurability(p)
                || !is_rep(n, length){
                continue;
            }

            let mut temp = n;
            let mut i = 0;
            loop{
                indices.insert(temp, (idx, i));
                temp = cyclic_move_unsafe(temp, length);
                i += 1;
                if temp == n{
                    break;
                }
            }

            basis.push((n, length));
            idx += 1;
        }

        if basis.len() == 0{
            return Err(Error::make_error_syntax(ErrorCode::InvalidConfiguration));
        }
        return Ok((basis, indices));
    }
}


#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_basis_nk(){
        let length = 4;
        let n = 0;
        let k = 0;

        let gen = BasisNK::new(EigenNumMomentum::new(n, k), length);
        let (base, indices) = gen.build().unwrap();
        assert_eq!(base, vec![NumMomentumState::new(SimpleState{rep : 0, length}, k).unwrap()]);
        for (&k, &v) in indices.iter(){
            assert!(k == 0 && v == (0, 0));
        }

        let n = 2;
        let gen = BasisNK::new(EigenNumMomentum::new(n, k), length);
        let (base, indices) = gen.build().unwrap();
        assert_eq!(base,
            vec![NumMomentumState::new(SimpleState{rep : 3, length}, k).unwrap(),
                NumMomentumState::new(SimpleState{rep : 5, length}, k).unwrap()]);
        for (&k, &v) in indices.iter(){
            assert!((k == 3 && v == (0, 0)) || (k == 5 && v == (1, 0)) || (k == 6 && v == (0, 3))
                    || (k == 9 && v == (0, 1)) || (k == 10 && v == (1, 1)) || (k == 12 && v == (0, 2)));
        }
    }

    #[test]
    fn test_basis_light_nk(){
        let length = 4;
        let n = 0;
        let k = 0;

        let gen = BasisNK::new(EigenNumMomentum::new(n, k), length);
        let (base, indices) = gen.build_light().unwrap();
        assert_eq!(base, vec![(0, 4)]);
        for (&k, &v) in indices.iter(){
            assert!(k == 0 && v == (0, 0));
        }

        let n = 2;
        let gen = BasisNK::new(EigenNumMomentum::new(n, k), length);
        let (base, indices) = gen.build_light().unwrap();
        assert_eq!(base, vec![(3, 4), (5, 4)]);
        for (&k, &v) in indices.iter(){
            assert!((k == 3 && v == (0, 0)) || (k == 5 && v == (1, 0)) || (k == 6 && v == (0, 3))
                    || (k == 9 && v == (0, 1)) || (k == 10 && v == (1, 1)) || (k == 12 && v == (0, 2)));
        }
    }
}
