
use crate::states::GeneralState;
use crate::prelude::*;


impl SimpleState{
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

#[derive(Copy, Clone, Debug, Hash, Eq, Ord, PartialEq, PartialOrd)]
pub struct EigenNumMomentum(usize, usize);
impl EigenValue for EigenNumMomentum {}
impl EigenNumMomentum{
    pub fn new(m : usize, k : usize) -> Self{
        EigenNumMomentum(m, k)
    }

    pub fn from_number(n : EigenNumber, k : usize) -> Self{
        EigenNumMomentum(n.total_number(), k)
    }

    pub fn total_number(&self) -> usize{
        self.0
    }

    pub fn wave_number(&self) -> usize{
        self.1
    }

    pub fn phase_factor(&self, length : usize) -> Complex64{
        Complex64::new(0f64, 2f64 * PI * (self.1 as f64) / (length as f64)).exp()
    }

    pub fn check_commensurability(&self, period : usize, length : usize) -> bool{
        (self.1 * period) % length == 0
    }
}
pub type NumMomentumState = EigenState<EigenNumMomentum>;

impl NumMomentumState{
    pub fn new(s : SimpleState, k : usize) -> Option<Self>{
        let (mut temp, period) = s.find_representative();
        let (rep, length) = (s.rep, s.length);
        let m = s.total_number();
        let eigen_v = EigenNumMomentum(m, k);

        if temp.rep != s.rep || !eigen_v.check_commensurability(period, length){
            return None;
        }

        let omega : Complex64 = eigen_v.phase_factor(length).inv();
        let mut coeff : Complex64 = Complex64::from( (period as f64).sqrt() / (length as f64));
        let mut result = NumMomentumState{
            state : GeneralState{
                states : vec![],
                coeffs : Array1::<Complex64>::zeros([period]),
            },
            index : RepWith::<EigenNumMomentum>::new(eigen_v, rep),
            length,
        };

        for i in 0..period{
            result.state.states.push(temp);
            result.state.coeffs[i] = coeff;

            temp.cyclic_move_mut();
            coeff *= omega;
        }


        return Some(result);
    }

    pub fn new_unsafe(s : SimpleState, k : usize) -> Self{
        let m = s.total_number();
        let p = s.period();
        let egn_nk = EigenNumMomentum(m, k);

        let omega = egn_nk.phase_factor(s.length).inv();
        let mut coeff = Complex64::from((p as f64).sqrt() / (s.length as f64));
        let mut result = NumMomentumState{
            state : GeneralState{
                states : vec![],
                coeffs : Array1::<Complex64>::zeros([p]),
            },
            index : RepWith::<EigenNumMomentum>::new(egn_nk, s.rep),
            length : s.length,
        };

        let mut temp = s.clone();
        for i in 0..p{
            result.state.states.push(temp);
            result.state.coeffs[i] = coeff;

            temp.cyclic_move_mut();
            coeff *= omega;
        }

        return result;
    }

    pub fn total_number(&self) -> usize{
        self.index.eigenvalue().0
    }

    pub fn wave_number(&self) -> usize{
        self.index.eigenvalue().1
    }

    pub fn phase_factor(&self) -> Complex64{
        Complex64::new(0f64, 2f64 * PI * (self.wave_number() as f64) / (self.length as f64)).exp()
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

    #[test]
    fn test_state() -> Result<(), Error>{
        let nkstate = NumMomentumState::new(SimpleState::new(5, 4), 0).unwrap();
        assert_eq!(nkstate.state.states[0], SimpleState::new(5, 4));
        assert_eq!(nkstate.state.states[1], SimpleState::new(10, 4));

        let c = Complex64::from(1f64 / 8f64.sqrt());
        aclose(nkstate.state.coeffs[[0]], c, 1e-10);
        aclose(nkstate.state.coeffs[[1]], c, 1e-10);

        assert_eq!(nkstate.index.eigenvalue().total_number(), 2);
        assert_eq!(nkstate.index.eigenvalue().wave_number(), 0);
        assert_eq!(nkstate.length, 4);

        let nkstate = NumMomentumState::new(SimpleState::new(10, 4), 0);
        assert_eq!(nkstate, None);

        let nkstate = NumMomentumState::new(SimpleState::new(5, 4), 1);
        assert_eq!(nkstate, None);
        Ok(())
    }
}
