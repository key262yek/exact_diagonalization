use crate::states::{SimpleState, number::{EigenNumber, NumberState}};
use fnv::FnvHashMap;
use num::integer::binomial;
use super::BasisGenerator;

pub type BasisN = BasisGenerator<EigenNumber>;
impl BasisN{
    pub fn new(v : EigenNumber, length : usize) -> Self{
        Self{
            length,
            value : Box::new(v),
        }
    }

    pub fn build(&self) -> (Vec<NumberState>, FnvHashMap<usize, usize>){
        let num = self.value.total_number();
        let length = self.length;

        let max_state = 1 << length;
        let mut basis : Vec<NumberState> = Vec::with_capacity(binomial(length, num));
        let mut indices : FnvHashMap<usize, usize> = FnvHashMap::default();
        let mut idx = 0;

        for n in 0..max_state {
            let state = SimpleState::new(n, length);
            if state.total_number() == num {
                basis.push(NumberState::new(state));
                indices.insert(n, idx);

                idx += 1;
            }
        }

        return (basis, indices);
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
        let (base, mut indices) = gen.build();
        assert_eq!(base, vec![NumberState::new(SimpleState{rep : 0, length})]);
        for (k, v) in indices.drain(){
            assert!(k == 0 && v == 0);
        }

        let rep = 2;
        let gen = BasisN::new(EigenNumber::new(rep), length);
        let (base, mut indices) = gen.build();
        assert_eq!(base,
            vec![NumberState::new(SimpleState{rep : 3, length}),
                 NumberState::new(SimpleState{rep : 5, length}),
                 NumberState::new(SimpleState{rep : 6, length}),
                 NumberState::new(SimpleState{rep : 9, length}),
                 NumberState::new(SimpleState{rep : 10, length}),
                 NumberState::new(SimpleState{rep : 12, length})]);
        for (k, v) in indices.drain(){
            assert!((k == 3 && v == 0) || (k == 5 && v == 1) || (k == 6 && v == 2)
                    || (k == 9 && v == 3) || (k == 10 && v == 4) || (k == 12 && v == 5));
        }
    }
}
