use std::hash::{Hash, Hasher};
use crate::prelude::*;

use std::cmp::Ordering;



pub mod bit_fns;
pub mod iterator;
pub mod number;
pub mod momentum;
pub mod symmetry;
pub mod representation;

pub trait State<T : EigenValue>{
    fn rep(&self) -> usize;
    fn length(&self) -> usize;
    fn value(&self) -> T;
    fn where_is(&self, num : usize) -> Option<usize>;
    // Position in linear combination of states with certain order

    fn coeff_of(&self, num : usize) -> Option<Complex64>;

    fn bit_sum(&self) -> usize{
        self.bit_iter().sum()
    }

    fn period(&self) -> usize{
        self.cycle_iter().count()
    }

    fn bit_iter(&self) -> BitIterator{
        match BitIterator::new(self.rep(), self.length()){
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        }
    }

    fn pair_iter(&self) -> PairIterator{
        match PairIterator::new(self.rep(), self.length()){
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        }
    }


    fn periodic_pair_iter(&self) -> PeriodicPairIterator{
        match PeriodicPairIterator::new(self.rep(), self.length()){
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        }
    }

    fn periodic_pair_enumerate(&self) -> PeriodicPairEnumerator{
        match PeriodicPairEnumerator::new(self.rep(), self.length()){
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        }
    }

    fn periodic_distanced_pair_iter(&self, dist : usize) -> PeriodicDistancedPairIterator{
        match PeriodicDistancedPairIterator::new(self.rep(), self.length(), dist){
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        }
    }

    fn periodic_distanced_pair_enumerate(&self, dist : usize) -> PeriodicDistancedPairEnumerator{
        match PeriodicDistancedPairEnumerator::new(self.rep(), self.length(), dist){
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        }
    }

    fn cycle_iter(&self) -> CycleIterator{
        match CycleIterator::new(self.rep(), self.length()){
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        }
    }
}

// =====================================================================================================
// =====================================================================================================


impl State<EmptyValue> for (usize, usize){
    fn rep(&self) -> usize{
        self.0
    }


    fn length(&self) -> usize{
        self.1
    }

    fn value(&self) -> EmptyValue{
        EmptyValue{}
    }

    fn where_is(&self, num : usize) -> Option<usize> {
        if self.0 == num{
            return Some(0);
        } else {
            None
        }
    }

    fn coeff_of(&self, num : usize) -> Option<Complex64> {
        if self.0 == num {
            return Some(Complex64::from(1.0));
        } else {
            None
        }
    }
}

// =====================================================================================================
// =====================================================================================================


#[derive(Copy, Clone, Debug, Eq, PartialEq, PartialOrd, Hash)]
pub struct SimpleState{
    pub rep : usize,
    pub length : usize,
}

impl SimpleState{
    pub fn new(rep : usize, length : usize) -> Self{
        if rep >= (1 << length) {
            panic!("{}", Error::make_error_syntax(ErrorCode::OverFlow));
        }

        Self{
            rep,
            length,
        }
    }


    pub fn period(&self) -> usize{
        return self.cycle_iter().count();
    }

    pub fn is_rep(&self) -> bool{
        for n in self.cycle_iter(){
            if n < self.rep {
                return false;
            }
        }
        return true;
    }

    pub fn find_representative(&self) -> (Self, usize){
        let mut period = 0;
        let mut min = self.rep;

        for n in self.cycle_iter(){
            period += 1;
            if n < min {
                min = n;
            }
        }
        return (Self{rep : min, length : self.length}, period);
    }

    pub fn way_to_representative(&self, rep : usize) -> Option<usize>{
        for (i, n) in self.cycle_iter().enumerate(){
            if n == rep {
                return Some(i);
            }
        }
        return None;
    }
}

impl State<EmptyValue> for SimpleState{
    fn rep(&self) -> usize{
        self.rep
    }


    fn length(&self) -> usize{
        self.length
    }

    fn value(&self) -> EmptyValue{
        EmptyValue{}
    }

    fn where_is(&self, num : usize) -> Option<usize> {
        if self.rep == num{
            return Some(0);
        } else {
            None
        }
    }

    fn coeff_of(&self, num : usize) -> Option<Complex64> {
        if self.rep == num {
            return Some(Complex64::from(1.0));
        } else {
            None
        }
    }
}


// =====================================================================================================
// =====================================================================================================



// =====================================================================================================
// =====================================================================================================

// pub trait AddInfo<W: ?Sized>{
//     type Info;
//     fn add_info(&self, info : Self::Info) -> W;
// }

// pub trait ConvertFrom<T>{
//     type Info;
//     fn from(value : T, info : Self::Info) -> Self;
// }


// =====================================================================================================
// =====================================================================================================


// impl AddInfo<EigenNumber> for EmptyValue{
//     type Info = usize;
//     fn add_info(&self, info : Self::Info) -> EigenNumber{
//         EigenNumber::new(info)
//     }
// }

// impl AddInfo<EigenNumMomentum> for EmptyValue{
//     type Info = (usize, usize);
//     fn add_info(&self, info : Self::Info) -> EigenNumMomentum{
//         EigenNumMomentum::new(info.0, info.1)
//     }
// }

// =====================================================================================================
// =====================================================================================================




// impl<T, W> AddInfo<Representation<W>> for Representation<T>
//     where T : AddInfo<W> + EigenValue,
//           W : EigenValue{
//     type Info = T::Info;
//     fn add_info(&self, info : Self::Info) -> Representation<W>{
//         Representation(self.0.add_info(info), self.1)
//     }
// }

// impl<T, W> ConvertFrom<Representation<W>> for Representation<T>
//     where T : ConvertFrom<W> + EigenValue,
//           W : EigenValue{
//     type Info = T::Info;
//     fn from(value : Representation<W>, info : Self::Info) -> Self{
//         Representation(T::from(value.0, info), value.1)
//     }
// }

// =====================================================================================================
// =====================================================================================================

#[derive(Clone, Debug)]
pub struct EigenState<T>
    where T : EigenValue{
    pub state : FnvHashMap<usize, (usize, Complex64)>,
    pub index : Representation<T>,
    pub length : usize,
}

impl<T> EigenState<T>
    where T : EigenValue{
    pub fn index(&self) -> Representation<T>{
        self.index
    }

    pub fn get_eigenvalue(&self) -> &T{
        self.index.get_eigenvalue()
    }

    pub fn get_mut_eigenvalue(&mut self) -> &mut T{
        self.index.get_mut_eigenvalue()
    }

    pub fn state(&self) -> &FnvHashMap<usize, (usize, Complex64)>{
        &self.state
    }
}

impl<T> PartialEq for EigenState<T>
    where T : EigenValue{
    fn eq(&self, other : &Self) -> bool{
        self.index == other.index
    }
}

impl<T> Eq for EigenState<T> where T : EigenValue + Copy + PartialEq{}

impl<T> PartialOrd for EigenState<T>
    where T : EigenValue + Ord{
    fn partial_cmp(&self, other : &Self) -> Option<Ordering>{
        Some(self.index.cmp(&other.index))
    }
}

impl<T> Ord for EigenState<T>
    where T : EigenValue + Ord{
    fn cmp(&self, other : &Self) -> Ordering{
        self.index.cmp(&other.index)
    }
}

impl<T> Hash for EigenState<T>
    where T : EigenValue{
    fn hash<H : Hasher>(&self, state : &mut H){
        self.index.hash(state);
    }
}

impl<T> State<T> for EigenState<T>
    where T : EigenValue{
    fn rep(&self) -> usize{
        self.index.get_rep()
    }

    fn length(&self) -> usize{
        self.length
    }

    fn value(&self) -> T{
        self.index.0
    }

    fn where_is(&self, num : usize) -> Option<usize>{
        self.state.get(&num).map(|x| x.0)
    }

    fn coeff_of(&self, num : usize) -> Option<Complex64> {
        self.state.get(&num).map(|x| x.1)
    }
}


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_period(){
        let s = SimpleState::new(0, 4);
        assert_eq!(s.period(), 1);

        let s = SimpleState::new(1, 4);
        assert_eq!(s.period(), 4);

        let s = SimpleState::new(3, 4);
        assert_eq!(s.period(), 4);

        let s = SimpleState::new(5, 4);
        assert_eq!(s.period(), 2);

        let s = SimpleState::new(7, 4);
        assert_eq!(s.period(), 4);

        let s = SimpleState::new(15, 4);
        assert_eq!(s.period(), 1);
    }

    #[test]
    fn test_is_rep(){
        let s = SimpleState::new(0, 4);
        assert!(s.is_rep());

        let s = SimpleState::new(1, 4);
        assert!(s.is_rep());

        let s = SimpleState::new(3, 4);
        assert!(s.is_rep());

        let s = SimpleState::new(5, 4);
        assert!(s.is_rep());

        let s = SimpleState::new(7, 4);
        assert!(s.is_rep());

        let s = SimpleState::new(15, 4);
        assert!(s.is_rep());

        let s = SimpleState::new(2, 4);
        assert!(!s.is_rep());

        let s = SimpleState::new(4, 4);
        assert!(!s.is_rep());

        let s = SimpleState::new(8, 4);
        assert!(!s.is_rep());
    }

    #[test]
    fn test_cycle_rep(){
        let s = SimpleState::new(0, 4);
        assert_eq!(s.find_representative(), (s, 1));

        let s = SimpleState::new(1, 4);
        assert_eq!(s.find_representative(), (s, 4));

        let s = SimpleState::new(3, 4);
        assert_eq!(s.find_representative(), (s, 4));

        let s = SimpleState::new(5, 4);
        assert_eq!(s.find_representative(), (s, 2));

        let s = SimpleState::new(7, 4);
        assert_eq!(s.find_representative(), (s, 4));

        let s = SimpleState::new(15, 4);
        assert_eq!(s.find_representative(), (s, 1));
    }

    #[test]
    fn test_way_to_rep(){
        let s = SimpleState::new(1, 4);
        assert_eq!(s.way_to_representative(1), Some(0));

        let s = SimpleState::new(2, 4);
        assert_eq!(s.way_to_representative(1), Some(1));

        let s = SimpleState::new(4, 4);
        assert_eq!(s.way_to_representative(1), Some(2));

        let s = SimpleState::new(8, 4);
        assert_eq!(s.way_to_representative(1), Some(3));

        let s = SimpleState::new(8, 4);
        assert_eq!(s.way_to_representative(3), None);
    }

    #[test]
    fn test_commensurability(){
        fn sum_commensurability(length : usize, period : usize){
            for k in 0..length{
                let omega : Complex64 = Complex64::new(0f64, -2f64 * PI * (k as f64) / ((length / period) as f64)).exp();
                let mut coeff : Complex64 = Complex64::from( 1.0);
                let mut sum : Complex64 = Complex64::from(0.0);
                let test : Complex64 = Complex64::from((length / period) as f64);

                for _ in 0..(length / period){
                    sum += coeff;
                    coeff *= omega;
                }

                if k * period % length == 0 {
                    assert!((sum - test).norm() < 1e-10);
                } else {
                    assert!(sum.norm() < 1e-10);
                }
            }
        }

        sum_commensurability(10, 2);
        sum_commensurability(10, 5);
        sum_commensurability(10, 1);
    }
}
