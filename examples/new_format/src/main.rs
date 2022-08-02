use std::collections::HashMap;
use num_complex::Complex64;

trait EigenValue {}

#[derive(Copy, Clone, Debug)]
struct NumberTrans(usize, usize);

impl EigenValue for NumberTrans {}

trait State<T : EigenValue>{
    fn eigenvalue(&self) -> T;
    fn rep(&self) -> usize;
    fn length(&self) -> usize;
    fn where_is(&self, other_rep : usize) -> Option<usize>;
}

struct TransState{
    rep : usize,
    length : usize,
    value : NumberTrans,
    superposition : HashMap<usize, (usize, Complex64)>
}

impl State<NumberTrans> for TransState{
    fn eigenvalue(&self) -> NumberTrans{
        self.value
    }

    fn rep(&self) -> usize{
        self.rep
    }

    fn length(&self) -> usize{
        self.length
    }

    fn where_is(&self, other_rep : usize) -> Option<usize>{
        self.superposition.get(&other_rep).map(|x| x.0)
    }
}



struct BasisGen;

impl BasisGen {
    fn build<S, T>(&self) -> (HashMap<T, Vec<S>>, HashMap<usize, (T, usize)>)
        where S : State<T>,
              T : EigenValue{
        unimplemented!();
    }
}

fn main(){
    let x = BasisGen;
    let (bases, indices) = x.build::<TransState, NumberTrans>();

    let xxz = PeriodicNearestXXZ::new(1f64, delta);
    let omega_k = basis[0].phase_factor();

    let n = basis.len();
    let mut hamiltonian : Array2<Complex64> = Array2::zeros((n, n));

    for (idx, state) in basis.iter().enumerate(){
        for (rep2, value) in xxz.apply_to(state){
            if let Some(idx2) = indices.get(&rep2){
                d =
                hamiltonian[[*idx2, idx]] += Complex64::from(value) * normal_f1 / normal_f2 * omega_k.powu(*d as u32);
            }
        }
    }
}

