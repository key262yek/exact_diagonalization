use crate::{prelude::*, states::bit_fns::{is_rep, sum_bit}};

pub mod number;
pub mod momentum;

#[derive(Clone, Debug)]
pub struct BasisGenerator<I : EigenValue>{
    pub length : usize,
    pub value : Box<I>,
}

pub type Basis = BasisGenerator<EmptyValue>;

// pub trait BasisGen<T : EigenValue>{
//     fn build(&self) -> (FnvHashMap<T, Vec<EigenState<T>>>, FnvHashMap<Representation<T>, usize>);
// }

// impl BasisGen<EigenNumber> for Basis{
//     fn build(&self) -> (FnvHashMap<EigenNumber, Vec<NumberState>>, FnvHashMap<RepNum, usize>){
//         let length = self.length;
//         let max_state = 1 << length;
//         let mut bases : FnvHashMap<EigenNumber, Vec<NumberState>> = FnvHashMap::default();
//         let mut indices : FnvHashMap<RepNum, usize> = FnvHashMap::default();

//        for n in 0..max_state {
//             let state = SimpleState::new(n, length);
//             let m = state.bit_sum();
//             let egn_v = EigenNumber::new(m);

//             match bases.get_mut(&egn_v){
//                 None => {
//                     let mut basis_n : Vec<NumberState> = Vec::with_capacity(binomial(length, m));
//                     basis_n.push(NumberState::new(&state));
//                     bases.insert(egn_v, basis_n);

//                     indices.insert(Representation::new(egn_v, n), 0);
//                 },
//                 Some(basis_n) => {
//                     let idx = basis_n.len();
//                     basis_n.push(NumberState::new(&state));
//                     indices.insert(Representation::new(egn_v, n), idx);
//                 },
//             }
//         }

//         return (bases, indices);
//     }
// }

// impl BasisGen<EigenNumMomentum> for Basis{
//     fn build(&self) -> (FnvHashMap<EigenNumMomentum, Vec<NumMomentumState>>, FnvHashMap<Representation<EigenNumMomentum>, usize>){

//     }
// }

impl Basis {
    pub fn new(l : usize) -> Self{
        Self{
            length : l,
            value : Box::new(EmptyValue::new())
        }
    }

    pub fn build_n(&self) -> (FnvHashMap<EigenNumber, Vec<NumberState>>, FnvHashMap<RepNum, (usize, usize)>){
        let length = self.length;
        let max_state = 1 << length;
        let mut bases : FnvHashMap<EigenNumber, Vec<NumberState>> = FnvHashMap::default();
        let mut indices : FnvHashMap<RepNum, (usize, usize)> = FnvHashMap::default();

       for n in 0..max_state {
            let state = SimpleState::new(n, length);
            let m = state.bit_sum();
            let egn_v = EigenNumber::new(m);

            match bases.get_mut(&egn_v){
                None => {
                    let mut basis_n : Vec<NumberState> = Vec::with_capacity(binomial(length, m));
                    basis_n.push(NumberState::new(&state));
                    bases.insert(egn_v, basis_n);

                    indices.insert(Representation::new(egn_v, n), (0, 0));
                },
                Some(basis_n) => {
                    let idx = basis_n.len();
                    basis_n.push(NumberState::new(&state));
                    indices.insert(Representation::new(egn_v, n), (idx, 0));
                },
            }
        }

        return (bases, indices);
    }

    pub fn build_light_n(&self) -> (FnvHashMap<EigenNumber, Vec<(usize, usize)>>, FnvHashMap<(EigenNumber, usize), (usize, usize)>){
        let length = self.length;
        let max_state = 1 << length;
        let mut bases : FnvHashMap<EigenNumber, Vec<(usize, usize)>> = FnvHashMap::default();
        let mut indices : FnvHashMap<(EigenNumber, usize), (usize, usize)> = FnvHashMap::default();

       for n in 0..max_state {
            let m = sum_bit(n);
            let egn_v = EigenNumber::new(m);

            match bases.get_mut(&egn_v){
                None => {
                    let mut basis_n : Vec<(usize, usize)> = Vec::with_capacity(binomial(length, m));
                    basis_n.push((n, length));
                    bases.insert(egn_v, basis_n);

                    indices.insert((egn_v, n), (0, 0));
                },
                Some(basis_n) => {
                    let idx = basis_n.len();
                    basis_n.push((n, length));
                    indices.insert((egn_v, n), (idx, 0));
                },
            }
        }

        return (bases, indices);
    }

    pub fn build_nk(&self) -> (FnvHashMap<EigenNumMomentum, Vec<NumMomentumState>>, FnvHashMap<RepNumMomentum, (usize, usize)>){

        let max_state = 1 << self.length;
        let mut bases : FnvHashMap<EigenNumMomentum, Vec<NumMomentumState>> = FnvHashMap::default();
        let mut indices : FnvHashMap<RepNumMomentum, (usize, usize)> = FnvHashMap::default();

        for n in 0..max_state {
            let state = SimpleState::new(n, self.length);
            if !state.is_rep(){
                continue;
            }

            let m = state.bit_sum();
            let period = state.period();
            for k in CommenIterator::new(period, self.length){
                let nkstate = NumMomentumState::new_unsafe(&state, &EigenNumMomentum(m, k));
                let egn_nk = EigenNumMomentum::new(m, k);

                match bases.get_mut(&egn_nk){
                    None => {
                        for (num, (i, _coeff)) in nkstate.state.iter(){
                            indices.insert(Representation(egn_nk, *num), (0, *i));
                        }

                        let mut basis_nk : Vec<NumMomentumState> = Vec::new();
                        basis_nk.push(nkstate);
                        bases.insert(egn_nk, basis_nk);

                    },
                    Some(basis_nk) => {
                        let idx = basis_nk.len();

                        for (num, (i, _coeff)) in nkstate.state.iter(){
                            indices.insert(Representation(egn_nk, *num), (idx, *i));
                        }
                        basis_nk.push(nkstate);
                    }
                };
            }
        }

        return (bases, indices);
    }

    pub fn build_light_nk(&self) -> (FnvHashMap<EigenNumMomentum, Vec<(usize, usize)>>, FnvHashMap<(EigenNumMomentum, usize), (usize, usize)>){

        let length = self.length;
        let max_state = 1 << length;
        let mut bases : FnvHashMap<EigenNumMomentum, Vec<(usize, usize)>> = FnvHashMap::default();
        let mut indices : FnvHashMap<(EigenNumMomentum, usize), (usize, usize)> = FnvHashMap::default();

        for n in 0..max_state {
            if !is_rep(n, length){
                continue;
            }

            let m = sum_bit(n);
            let period = period_unsafe(n, length);
            for k in CommenIterator::new(period, length){
                let egn_nk = EigenNumMomentum::new(m, k);

                match bases.get_mut(&egn_nk){
                    None => {
                        let mut temp = n;
                        let mut i = 0;
                        loop{
                            indices.insert((egn_nk, temp), (0, i));
                            temp = cyclic_move_unsafe(temp, length);
                            i += 1;
                            if temp == n{
                                break;
                            }
                        }

                        bases.insert(egn_nk, vec![(n, length)]);

                    },
                    Some(basis_nk) => {
                        let idx = basis_nk.len();

                        let mut temp = n;
                        let mut i = 0;
                        loop{
                            indices.insert((egn_nk, temp), (idx, i));
                            temp = cyclic_move_unsafe(temp, length);
                            i += 1;
                            if temp == n{
                                break;
                            }
                        }

                        basis_nk.push((n, length));
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
                        vec![NumberState::new(&SimpleState::new(0, length))]);
                },
                1 => {
                    assert_eq!(*basis_n,
                        vec![NumberState::new(&SimpleState::new(1, length)),
                             NumberState::new(&SimpleState::new(2, length))]);
                },
                2 => {
                    assert_eq!(*basis_n,
                        vec![NumberState::new(&SimpleState::new(3, length))]);
                },
                _ => unreachable!()
            }
        }

        for (Representation(egn_v, n), idx) in indices.iter(){
            match n{
                0 => {
                    assert_eq!(egn_v.total_number(), 0);
                    assert_eq!(idx.0, 0);
                },
                1 => {
                    assert_eq!(egn_v.total_number(), 1);
                    assert_eq!(idx.0, 0);
                },
                2 => {
                    assert_eq!(egn_v.total_number(), 1);
                    assert_eq!(idx.0, 1);
                },
                3 => {
                    assert_eq!(egn_v.total_number(), 2);
                    assert_eq!(idx.0, 0);
                },
                _ => unreachable!()
            }
        }
    }

    #[test]
    fn test_whole_basis_light_n(){
        let length = 2;
        let basis_gen  = Basis::new(length);
        let (base, indices) = basis_gen.build_light_n();

        for (egn_v, basis_n) in base.iter(){
            match egn_v.total_number(){
                0 => {
                    assert_eq!(*basis_n, vec![(0, 2)]);
                },
                1 => {
                    assert_eq!(*basis_n, vec![(1, 2), (2, 2)]);
                },
                2 => {
                    assert_eq!(*basis_n, vec![(3, 2)]);
                },
                _ => unreachable!()
            }
        }

        for ((egn_v, n), idx) in indices.iter(){
            match n{
                0 => {
                    assert_eq!(egn_v.total_number(), 0);
                    assert_eq!(idx.0, 0);
                },
                1 => {
                    assert_eq!(egn_v.total_number(), 1);
                    assert_eq!(idx.0, 0);
                },
                2 => {
                    assert_eq!(egn_v.total_number(), 1);
                    assert_eq!(idx.0, 1);
                },
                3 => {
                    assert_eq!(egn_v.total_number(), 2);
                    assert_eq!(idx.0, 0);
                },
                _ => unreachable!()
            }
        }
    }

    #[test]
    fn test_whole_basis_nk(){
        let length = 4;
        let egn_v = EigenNumMomentum(2, 0);
        let basis_gen  = Basis::new(length);
        let (base, indices) = basis_gen.build_nk();

        assert_eq!(base.get(&EigenNumMomentum::new(2, 0)),
            Some(&vec![NumMomentumState::new(&(3, length), &egn_v).unwrap(),
                NumMomentumState::new(&(5, length), &egn_v).unwrap()]));

        assert_eq!(indices.get(&Representation(egn_v, 3)), Some(&(0, 0)));
        assert_eq!(indices.get(&Representation(egn_v, 5)), Some(&(1, 0)));
        assert_eq!(indices.get(&Representation(egn_v, 6)), Some(&(0, 3)));
        assert_eq!(indices.get(&Representation(egn_v, 9)), Some(&(0, 1)));
        assert_eq!(indices.get(&Representation(egn_v, 10)), Some(&(1, 1)));
        assert_eq!(indices.get(&Representation(egn_v, 12)), Some(&(0, 2)));
        assert_eq!(indices.get(&Representation(egn_v, 0)), None);
    }

    #[test]
    fn test_whole_basis_light_nk(){
        let length = 4;
        let basis_gen  = Basis::new(length);
        let (base, indices) = basis_gen.build_light_nk();

        assert_eq!(base.get(&EigenNumMomentum::new(2, 0)), Some(&vec![(3, 4), (5, 4)]));

        assert_eq!(indices.get(&(EigenNumMomentum::new(2, 0), 3)), Some(&(0, 0)));
        assert_eq!(indices.get(&(EigenNumMomentum::new(2, 0), 5)), Some(&(1, 0)));
        assert_eq!(indices.get(&(EigenNumMomentum::new(2, 0), 6)), Some(&(0, 3)));
        assert_eq!(indices.get(&(EigenNumMomentum::new(2, 0), 9)), Some(&(0, 1)));
        assert_eq!(indices.get(&(EigenNumMomentum::new(2, 0), 10)), Some(&(1, 1)));
        assert_eq!(indices.get(&(EigenNumMomentum::new(2, 0), 12)), Some(&(0, 2)));
        assert_eq!(indices.get(&(EigenNumMomentum::new(2, 0), 0)), None);
    }
}
