use crate::prelude::*;
use super::{GeneralState, AddInfo};


#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct EigenNumber(usize);
impl EigenValue for EigenNumber {}
impl AddInfo<EigenNumMomentum> for EigenNumber{
    fn add_info(&self, info : usize) -> EigenNumMomentum{
        EigenNumMomentum::new(self.0, info)
    }
}

impl EigenNumber{
    pub fn new(n : usize) -> Self{
        Self(n)
    }

    pub fn total_number(&self) -> usize{
        self.0
    }
}

pub type NumberState = EigenState<EigenNumber>;

impl NumberState{
    pub fn new(s : SimpleState) -> Self{
        let m = s.total_number();
        Self{
            state : GeneralState{
                states : vec![s],
                coeffs : arr1(&[Complex64::new(1.0, 0.0)]),
            },
            index : RepWith::<EigenNumber>::new(EigenNumber(m), s.rep),
            length : s.length,
        }
    }

    pub fn from_rep(rep : usize, length : usize) -> Self{
        let s = SimpleState::new(rep, length);
        NumberState::new(s)
    }

    pub fn total_number(&self) -> usize{
        self.index.eigenvalue().0
    }
}

impl SimpleState{
    pub fn total_number(&self) -> usize{
        self.bit_iter().sum()
    }
}



#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_total_sum(){
        let state = SimpleState::new(10, 5);
        assert_eq!(state.total_number(), 2);

        let state = SimpleState::new(11, 5);
        assert_eq!(state.total_number(), 3);

        let state = SimpleState::new(15, 5);
        assert_eq!(state.total_number(), 4);
    }
}
