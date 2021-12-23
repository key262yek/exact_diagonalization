use crate::states::iterator::{BitIterator, PairIterator, PeriodicPairIterator, PeriodicPairEnumerator, CycleIterator};
use ndarray::Array1;
use num_complex::Complex;
use std::cmp::Ordering;
use std::hash::{Hash, Hasher};
use std::io::Empty;
use crate::error::{Error, ErrorCode};

use self::number::EigenNumber;

pub mod bit_fns;
pub mod iterator;
pub mod number;
pub mod momentum;


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
}

impl SimpleState{
    pub fn bit_iter(&self) -> BitIterator{
        match BitIterator::new(self.rep, self.length){
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        }
    }

    pub fn pair_iter(&self) -> PairIterator{
        match PairIterator::new(self.rep, self.length){
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        }
    }


    pub fn periodic_pair_iter(&self) -> PeriodicPairIterator{
        match PeriodicPairIterator::new(self.rep, self.length){
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        }
    }

    pub fn periodic_pair_enumerate(&self) -> PeriodicPairEnumerator{
        match PeriodicPairEnumerator::new(self.rep, self.length){
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        }
    }


    pub fn cycle_iter(&self) -> CycleIterator{
        match CycleIterator::new(self.rep, self.length){
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        }
    }
}

// =====================================================================================================
// =====================================================================================================

#[derive(Clone, Debug)]
pub struct GeneralState{
    pub states : Vec<SimpleState>,
    pub coeffs : Array1<Complex<f64>>,
}

// =====================================================================================================
// =====================================================================================================

pub trait EigenValue{}

// =====================================================================================================
// =====================================================================================================

pub trait AddInfo<W>{
    fn add_info(&self, info : usize) -> W;
}

// =====================================================================================================
// =====================================================================================================

pub struct EmptyValue();
impl EigenValue for EmptyValue{}
impl EmptyValue{
    pub fn new() -> Self{
        Self()
    }
}

impl AddInfo<EigenNumber> for EmptyValue{
    fn add_info(&self, info : usize) -> EigenNumber{
        EigenNumber::new(info)
    }
}

// =====================================================================================================
// =====================================================================================================

pub trait Representation<T : EigenValue>{
    fn value(&self) -> T;
    fn num(&self) -> usize;
}
impl Representation<EmptyValue> for usize {
    fn value(&self) -> EmptyValue{
        EmptyValue()
    }

    fn num(&self) -> usize{
        *self
    }
}

pub struct RepWith<T : EigenValue>(T, usize);
impl<T : EigenValue + Copy> Representation<T> for RepWith<T>{
    fn value(&self) -> T{
        self.0
    }
    fn num(&self) -> usize{
        self.1
    }
}

impl<T, W> AddInfo<RepWith<W>> for RepWith<T>
    where T : AddInfo<W> + EigenValue,
          W : EigenValue{
    fn add_info(&self, info : usize) -> RepWith<W>{
        RepWith(self.0.add_info(info), self.1)
    }
}

// =====================================================================================================
// =====================================================================================================

#[derive(Clone, Debug)]
pub struct EigenState<T>
    where T : EigenValue + Clone {
    pub state : GeneralState,
    pub rep : usize,
    pub length : usize,
    pub value : Box<T>,
}

impl<T> EigenState<T>
    where T : EigenValue + Clone{
    pub fn rep(&self) -> usize{
        self.rep
    }

    pub fn length(&self) -> usize{
        self.length
    }
}

impl<T> PartialEq for EigenState<T>
    where T : EigenValue + Clone + PartialEq{
    fn eq(&self, other : &Self) -> bool{
        self.value == other.value
    }
}

impl<T> Eq for EigenState<T> where T : EigenValue + Clone + PartialEq{}

impl<T> PartialOrd for EigenState<T>
    where T : EigenValue + Clone + Ord{
    fn partial_cmp(&self, other : &Self) -> Option<Ordering>{
        Some(self.value.cmp(&other.value))
    }
}

impl<T> Ord for EigenState<T>
    where T : EigenValue + Clone + Ord {
    fn cmp(&self, other : &Self) -> Ordering{
        self.value.cmp(&other.value)
    }
}

impl<T> Hash for EigenState<T>
    where T : EigenValue + Clone + Hash {
    fn hash<H : Hasher>(&self, state : &mut H){
        self.value.hash(state);
    }
}

impl<T> EigenState<T>
    where T : EigenValue + Clone{
    pub fn bit_iter(&self) -> BitIterator{
        match BitIterator::new(self.rep, self.length){
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        }
    }

    pub fn pair_iter(&self) -> PairIterator{
        match PairIterator::new(self.rep, self.length){
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        }
    }


    pub fn periodic_pair_iter(&self) -> PeriodicPairIterator{
        match PeriodicPairIterator::new(self.rep, self.length){
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        }
    }

    pub fn periodic_pair_enumerate(&self) -> PeriodicPairEnumerator{
        match PeriodicPairEnumerator::new(self.rep, self.length){
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        }
    }


    pub fn cycle_iter(&self) -> CycleIterator{
        match CycleIterator::new(self.rep, self.length){
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        }
    }
}
