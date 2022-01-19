use crate::prelude::*;

pub type NumberState = EigenState<EigenNumber>;

impl NumberState{
    pub fn new<S, T>(s : &S) -> Self
        where S : State<T>, T : EigenValue{
        let m = s.bit_sum();
        let mut state : FnvHashMap<usize, (usize, Complex64)> = FnvHashMap::default();
        state.insert(s.rep(), (0, Complex64::from(1f64)));
        Self{
            state : state,
            index : Representation(EigenNumber(m), s.rep()),
            length : s.length(),
        }
    }

    pub fn from_rep(rep : usize, length : usize) -> Self{
        NumberState::new(&(rep, length))
    }
}





#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_total_sum(){
        let state = SimpleState::new(10, 5);
        assert_eq!(state.bit_sum(), 2);

        let state = SimpleState::new(11, 5);
        assert_eq!(state.bit_sum(), 3);

        let state = SimpleState::new(15, 5);
        assert_eq!(state.bit_sum(), 4);
    }
}
