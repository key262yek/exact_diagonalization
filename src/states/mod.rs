use crate::prelude::*;

use std::cmp::Ordering;
use std::hash::{Hash, Hasher};


pub mod bit_fns;
pub mod iterator;
pub mod number;
pub mod momentum;

pub trait State<T : EigenValue>{
    fn rep(&self) -> usize;
    fn length(&self) -> usize;
    fn value(&self) -> T;

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
}

// =====================================================================================================
// =====================================================================================================

#[derive(Clone, Debug)]
pub struct GeneralState{
    pub states : Vec<SimpleState>,
    pub coeffs : Array1<Complex64>,
}

// =====================================================================================================
// =====================================================================================================

pub trait EigenValue : Copy + Eq + PartialEq + Hash{}

impl EigenValue for usize {}

// =====================================================================================================
// =====================================================================================================

pub trait AddInfo<W: ?Sized>{
    type Info;
    fn add_info(&self, info : Self::Info) -> W;
}

pub trait ConvertFrom<T>{
    type Info;
    fn from(value : T, info : Self::Info) -> Self;
}


// =====================================================================================================
// =====================================================================================================

#[derive(Copy, Clone, Hash, Debug, Eq, PartialEq)]
pub struct EmptyValue();
impl EigenValue for EmptyValue{}
impl EmptyValue{
    pub fn new() -> Self{
        Self()
    }
}

impl AddInfo<EigenNumber> for EmptyValue{
    type Info = usize;
    fn add_info(&self, info : Self::Info) -> EigenNumber{
        EigenNumber::new(info)
    }
}

impl AddInfo<EigenNumMomentum> for EmptyValue{
    type Info = (usize, usize);
    fn add_info(&self, info : Self::Info) -> EigenNumMomentum{
        EigenNumMomentum::new(info.0, info.1)
    }
}

// =====================================================================================================
// =====================================================================================================

pub trait Representation<T : EigenValue>{
    fn eigenvalue(&self) -> T;
    fn rep(&self) -> usize;
}

impl Representation<EmptyValue> for usize {
    fn eigenvalue(&self) -> EmptyValue{
        EmptyValue()
    }

    fn rep(&self) -> usize{
        *self
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct RepWith<T : EigenValue>(pub T, pub usize);

impl<T> RepWith<T>
    where T : EigenValue{
    pub fn new(egn_v : T, rep : usize) -> Self{
        Self(egn_v, rep)
    }
}

impl<T : EigenValue> Representation<T> for RepWith<T>{
    fn eigenvalue(&self) -> T{
        self.0
    }
    fn rep(&self) -> usize{
        self.1
    }
}

impl<T, W> AddInfo<RepWith<W>> for RepWith<T>
    where T : AddInfo<W> + EigenValue,
          W : EigenValue{
    type Info = T::Info;
    fn add_info(&self, info : Self::Info) -> RepWith<W>{
        RepWith(self.0.add_info(info), self.1)
    }
}

impl<T, W> ConvertFrom<RepWith<W>> for RepWith<T>
    where T : ConvertFrom<W> + EigenValue,
          W : EigenValue{
    type Info = T::Info;
    fn from(value : RepWith<W>, info : Self::Info) -> Self{
        RepWith(T::from(value.0, info), value.1)
    }
}

// =====================================================================================================
// =====================================================================================================

#[derive(Clone, Debug)]
pub struct EigenState<T>
    where T : EigenValue{
    pub state : GeneralState,
    pub index : RepWith<T>,
    pub length : usize,
}

impl<T> EigenState<T>
    where T : EigenValue{
    pub fn index(&self) -> RepWith<T>{
        self.index
    }

    pub fn eigenvalue(&self) -> T{
        self.index.eigenvalue()
    }

    pub fn states(&self) -> &Vec<SimpleState>{
        &self.state.states
    }

    pub fn coeffs(&self) -> &Array1<Complex64>{
        &self.state.coeffs
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
        self.index.rep()
    }

    fn length(&self) -> usize{
        self.length
    }

    fn value(&self) -> T{
        self.index.0
    }
}
