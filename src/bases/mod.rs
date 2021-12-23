use num::integer::binomial;
use crate::states::{EigenValue, EmptyValue, SimpleState, iterator::CommenIterator, momentum::{EigenNumMomentum, NumMomentumState}, number::{EigenNumber, NumberState}};
use fnv::FnvHashMap;

pub mod number;
pub mod momentum;

#[derive(Clone, Debug)]
pub struct BasisGenerator<I : EigenValue>{
    pub length : usize,
    pub value : Box<I>,
}

pub type Basis = BasisGenerator<EmptyValue>;

impl Basis {
    pub fn new(l : usize) -> Self{
        Self{
            length : l,
            value : Box::new(EmptyValue::new())
        }
    }

    pub fn build_n(&self) -> (FnvHashMap<EigenNumber, Vec<NumberState>>, FnvHashMap<(EigenNumber, usize), usize>){
        let length = self.length;
        let max_state = 1 << length;
        let mut magnet_sets : FnvHashMap<EigenNumber, Vec<NumberState>> = FnvHashMap::default();
        let mut indices : FnvHashMap<(EigenNumber, usize), usize> = FnvHashMap::default();

       for n in 0..max_state {
            let state = SimpleState::new(n, length);
            let m = state.total_number();
            let egn_v = EigenNumber::new(m);

            match magnet_sets.get_mut(&egn_v){
                None => {
                    let mut basis_n : Vec<NumberState> = Vec::with_capacity(binomial(length, m));
                    basis_n.push(NumberState::new(state));
                    magnet_sets.insert(egn_v, basis_n);

                    indices.insert((egn_v, n), 0);
                },
                Some(basis_n) => {
                    let idx = basis_n.len();
                    basis_n.push(NumberState::new(state));
                    indices.insert((egn_v, n), idx);
                },
            }
        }

        return (magnet_sets, indices);
    }

    pub fn build_nk(&self) -> (FnvHashMap<EigenNumMomentum, Vec<NumMomentumState>>, FnvHashMap<(EigenNumMomentum, usize), (usize, usize)>){

        let max_state = 1 << self.length;
        let mut bases : FnvHashMap<EigenNumMomentum, Vec<NumMomentumState>> = FnvHashMap::default();
        let mut indices : FnvHashMap<(EigenNumMomentum, usize), (usize, usize)> = FnvHashMap::default();

        for n in 0..max_state {
            let state = SimpleState::new(n, self.length);
            if !state.is_rep(){
                continue;
            }

            let m = state.total_number();
            let period = state.period();
            for k in CommenIterator::new(period, self.length){
                let nkstate = NumMomentumState::new_unsafe(state, k);
                let egn_nk = EigenNumMomentum::new(m, k);

                match bases.get_mut(&egn_nk){
                    None => {
                        for (i, state) in nkstate.state.states.iter().enumerate(){
                            let num = state.rep;
                            indices.insert((egn_nk, num), (0, (period - i) % period));
                        }

                        let mut basis_nk : Vec<NumMomentumState> = Vec::new();
                        basis_nk.push(nkstate);
                        bases.insert(egn_nk, basis_nk);

                    },
                    Some(basis_nk) => {
                        let idx = basis_nk.len();

                        for (i, state) in nkstate.state.states.iter().enumerate(){
                            let num = state.rep;
                            indices.insert((egn_nk, num), (idx, (period - i) % period));
                        }
                        basis_nk.push(nkstate);
                    }
                };
            }
        }

        return (bases, indices);
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_whole_basis_n(){
        let length = 2;
        let basis_gen  = Basis::new(length);
        let (base, indices) = basis_gen.build_n();

        for (egn_v, basis_n) in base.iter(){
            match egn_v.total_number(){
                0 => {
                    assert_eq!(*basis_n,
                        vec![NumberState::new(SimpleState::new(0, length))]);
                },
                1 => {
                    assert_eq!(*basis_n,
                        vec![NumberState::new(SimpleState::new(1, length)),
                             NumberState::new(SimpleState::new(2, length))]);
                },
                2 => {
                    assert_eq!(*basis_n,
                        vec![NumberState::new(SimpleState::new(3, length))]);
                },
                _ => unreachable!()
            }
        }

        for ((egn_v, n), idx) in indices.iter(){
            match n{
                0 => {
                    assert_eq!(egn_v.total_number(), 0);
                    assert_eq!(*idx, 0);
                },
                1 => {
                    assert_eq!(egn_v.total_number(), 1);
                    assert_eq!(*idx, 0);
                },
                2 => {
                    assert_eq!(egn_v.total_number(), 1);
                    assert_eq!(*idx, 1);
                },
                3 => {
                    assert_eq!(egn_v.total_number(), 2);
                    assert_eq!(*idx, 0);
                },
                _ => unreachable!()
            }
        }
    }
}
