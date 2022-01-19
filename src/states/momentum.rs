
use crate::prelude::*;

use super::representation::FindRepresentation;

pub type NumMomentumState = EigenState<EigenNumMomentum>;

impl NumMomentumState{
    pub fn new<T>(s : &dyn State<T>, eigen_v : &EigenNumMomentum) -> Option<Self>
        where T : EigenValue + LowerThan<EigenNumMomentum>{
        // Return state only if s is representative state of eigen_v

        if !eigen_v.is_rep(s){
            return None;
        };
        let rep = s.rep();
        let length = s.length();
        let period = s.period();

        let omega : Complex64 = eigen_v.phase_factor(length).inv();
        let mut coeff : Complex64 = Complex64::from( (period as f64).sqrt() / (length as f64));

        let mut result = NumMomentumState{
            state : FnvHashMap::<usize, (usize, Complex64)>::default(),
            index : Representation(*eigen_v, rep),
            length,
        };

        for (idx, n) in (rep, length).cycle_iter().enumerate(){
            result.state.insert(n, (idx, coeff));
            coeff *= omega;
        }

        return Some(result);
    }

    pub fn new_unsafe<T>(s : &dyn State<T>, eigen_v : &EigenNumMomentum) -> Self
        where T : EigenValue + LowerThan<EigenNumMomentum>{
        // Use when s is representative state of EigenNumMomentum in certain.

        let period = s.period();
        let length = s.length();
        let rep = s.rep();


        let omega = eigen_v.phase_factor(length).inv();
        let mut coeff = Complex64::from((period as f64).sqrt() / (length as f64));
        let mut result = NumMomentumState{
            state : FnvHashMap::<usize, (usize, Complex64)>::default(),
            index : Representation(*eigen_v, rep),
            length : length,
        };

        for (idx, n) in (rep, length).cycle_iter().enumerate(){
            result.state.insert(n, (idx, coeff));
            coeff *= omega;
        }

        return result;
    }

    pub fn total_number(&self) -> usize{
        self.index.get_eigenvalue().total_number()
    }

    pub fn wave_number(&self) -> usize{
        self.index.get_eigenvalue().wave_number()
    }

    pub fn phase_factor(&self) -> Complex64{
        self.index.get_eigenvalue().phase_factor(self.length())
        // Complex64::new(0f64, 2f64 * PI * (self.wave_number() as f64) / (self.length as f64)).exp()
    }

    pub fn normalize_factor(&self) -> Complex64{
        let p = period_unsafe(self.rep(), self.length);
        Complex64::from((p as f64).sqrt() / (self.length as f64))
    }
}


#[cfg(test)]
mod test {
    use ndarray_linalg::{aclose};

    use super::*;

    #[test]
    fn test_state() -> Result<(), Error>{
        let nkstate = NumMomentumState::new(&SimpleState::new(5, 4),
                                                             &EigenNumMomentum(2, 0)).unwrap();
        assert_eq!(nkstate.where_is(5), Some(0));
        assert_eq!(nkstate.where_is(10), Some(1));

        let c = Complex64::from(1f64 / 8f64.sqrt());
        aclose(nkstate.coeff_of(5).unwrap(), c, 1e-10);
        aclose(nkstate.coeff_of(10).unwrap(), c, 1e-10);

        assert_eq!(nkstate.index.get_eigenvalue().total_number(), 2);
        assert_eq!(nkstate.index.get_eigenvalue().wave_number(), 0);
        assert_eq!(nkstate.length, 4);

        let nkstate = NumMomentumState::new(&SimpleState::new(10, 4),
                                                                   &EigenNumMomentum(2, 0));
        assert_eq!(nkstate, None);

        let nkstate = NumMomentumState::new(&SimpleState::new(5, 4),
                                                                         &EigenNumMomentum(2, 1));
        assert_eq!(nkstate, None);
        Ok(())
    }
}
