use crate::prelude::*;
use std::hash::Hash;


pub trait EigenValue : Copy + Eq + PartialEq + Hash{}

pub trait LowerThan<T : EigenValue> : EigenValue {
    fn check_extensible(&self, other : &T) -> bool;
}
pub trait HigherThan<T : EigenValue> : EigenValue {
    fn check_projectible(&self, other : &T) -> bool;
}


impl<T, W> HigherThan<T> for W
    where T : LowerThan<W> + EigenValue,
          W : EigenValue {
    fn check_projectible(&self, other : &T) -> bool {
        other.check_extensible(self)
    }
}

impl EigenValue for usize {}
impl<T> LowerThan<T> for usize where T : EigenValue{
    fn check_extensible(&self, _other : &T) -> bool {
        true
    }
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

impl<T> LowerThan<T> for EmptyValue where T : EigenValue{
    fn check_extensible(&self, _other : &T) -> bool {
        true
    }
}



// =====================================================================================================
// =====================================================================================================


pub trait NumberConservation : EigenValue{
    fn total_number(&self) -> usize;
}


#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct EigenNumber(pub usize);

impl EigenNumber{
    pub fn new(n : usize) -> Self{
        Self(n)
    }
}

impl EigenValue for EigenNumber {}

impl NumberConservation for  EigenNumber{
    fn total_number(&self) -> usize{
        self.0
    }
}

impl<T> LowerThan<T> for EigenNumber where T : EigenValue + NumberConservation {
    fn check_extensible(&self, other : &T) -> bool {
        self.total_number() == other.total_number()
    }
}

// =====================================================================================================
// =====================================================================================================


pub trait TranslationalSymmetry : EigenValue{
    fn wave_number(&self) -> usize;

    fn phase_factor(&self, length : usize) -> Complex64;

    fn check_commensurability(&self, period : usize, length : usize) -> bool;
}


#[derive(Copy, Clone, Debug, Hash, Eq, Ord, PartialEq, PartialOrd)]
pub struct EigenNumMomentum(pub usize,pub usize);
impl EigenValue for EigenNumMomentum {}

impl EigenNumMomentum{
    pub fn new(m : usize, k : usize) -> Self{
        EigenNumMomentum(m, k)
    }

    pub fn from_number(n : EigenNumber, k : usize) -> Self{
        EigenNumMomentum(n.total_number(), k)
    }
}

impl NumberConservation for EigenNumMomentum{
    fn total_number(&self) -> usize {
        self.0
    }
}

impl TranslationalSymmetry for EigenNumMomentum{
    fn wave_number(&self) -> usize{
        self.1
    }

    fn phase_factor(&self, length : usize) -> Complex64{
        Complex64::new(0f64, 2f64 * PI * (self.1 as f64) / (length as f64)).exp()
    }

    fn check_commensurability(&self, period : usize, length : usize) -> bool{
        (self.1 * period) % length == 0
    }
}

impl<T> LowerThan<T> for EigenNumMomentum where T : EigenValue + NumberConservation + TranslationalSymmetry {
    fn check_extensible(&self, other : &T) -> bool {
        (self.total_number() == other.total_number())
        && (self.wave_number() == other.wave_number())
    }
}

