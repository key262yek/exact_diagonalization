use crate::prelude::*;


#[derive(Clone, Copy, Debug, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct Representation<T : EigenValue>(pub T, pub usize);

pub type RepNum = Representation<EigenNumber>;
pub type RepNumMomentum = Representation<EigenNumMomentum>;

impl<T> Representation<T>
    where T : EigenValue{
    pub fn new(egn_v : T, rep : usize) -> Self{
        Self(egn_v, rep)
    }
}

impl<T : EigenValue> Representation<T>{
    pub fn get_eigenvalue(&self) -> &T{
        &self.0
    }

    pub fn get_rep(&self) -> usize{
        self.1
    }

    pub fn get_mut_eigenvalue(&mut self) -> &mut T{
        &mut self.0
    }

    pub fn get_mut_rep(&mut self) -> &mut usize{
        &mut self.1
    }
}

pub trait FindRepresentation<T : LowerThan<Self>> : EigenValue{
    fn find_rep(&self, state : &dyn State<T>) -> Option<usize>;
    fn is_rep(&self, state : &dyn State<T>) -> bool;
}

impl<T> FindRepresentation<T> for EmptyValue
    where T : EigenValue + LowerThan<EmptyValue>{
    fn find_rep(&self, state : &dyn State<T>) -> Option<usize> {
        Some(state.rep())
    }

    fn is_rep(&self, _state : &dyn State<T>) -> bool {
        true
    }
}

impl<T> FindRepresentation<T> for EigenNumber
    where T : EigenValue + LowerThan<EigenNumber>{
    fn find_rep(&self, state : &dyn State<T>) -> Option<usize> {
        Some(state.rep())
    }

    fn is_rep(&self, _state : &dyn State<T>) -> bool {
        true
    }
}

impl<T> FindRepresentation<T> for EigenNumMomentum
    where T : EigenValue + LowerThan<EigenNumMomentum>{
    fn find_rep(&self, state : &dyn State<T>) -> Option<usize> {
        let mut period = 0;
        let mut min = state.rep();

        for n in state.cycle_iter(){
            period += 1;
            if n < min {
                min = n;
            }
        }

        if self.check_commensurability(period, state.length()){
            Some(min)
        } else {
            None
        }
    }

    fn is_rep(&self, state : &dyn State<T>) -> bool {
        let mut period = 0;
        let rep = state.rep();

        for n in state.cycle_iter(){
            period += 1;
            if n < rep {
                return false;
            }
        }

        self.check_commensurability(period, state.length())
    }
}

// impl<T : EigenValue> FindRepresentation<T> for EigenState<T>{
//     fn find_rep(&self) -> Self {
//         let eig_v = self.index.get_eigenvalue();
//         let mut temp = self.clone();
//         *temp.index.get_mut_rep() = eig_v.find_rep_from(self);
//         return temp
//     }

//     fn is_rep(&self) -> bool {
//         self.index.get_eigenvalue().check_rep(self)
//     }

//     fn step_from(&self, rep : &Self) -> Option<usize> {
//         self.index.get_eigenvalue().step_from_rep(self, rep)
//     }

//     fn step_from_rep(&self) -> usize {
//         let rep = self.find_rep();
//         self.step_from(&rep).unwrap()
//     }
// }
