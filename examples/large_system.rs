use ndarray_linalg::{Eigh, UPLO};
use exact_diagonalization::prelude::*;

fn hamiltonian_with(l : usize, m : usize, k : usize) -> Result<(Array1<f64>, Array2<Complex64>), Error>{
    let basis_gen = BasisNK::new(EigenNumMomentum::new(m, k), l);
    let (basis, indices) = basis_gen.build();

    let delta = 2f64;
    let xxz = PeriodicNearestXXZ::new(delta);
    let omega_k = basis[0].phase_factor();

    let n = basis.len();
    let mut hamiltonian : Array2<Complex64> = Array2::zeros((n, n));

    for (idx, state) in basis.iter().enumerate(){
        let normal_f1 = state.normalize_factor();
        for (rep2, value) in xxz.apply_to(state){
            if let Some((idx2, d)) = indices.get(&rep2){
                let normal_f2 = basis[*idx2].normalize_factor();
                hamiltonian[[*idx2, idx]] += Complex64::from(value) * normal_f1 / normal_f2 * omega_k.powu(*d as u32);
            }
        }
    }

    return hamiltonian.eigh(UPLO::Lower).map_err(|c| Error::make_error_msg(format_args!("{}", c).to_string()));
}


fn main() -> Result<(), Error>{
    let length = 13;
    let m = 7;
    let k = 0;

    let (values, vectors) = hamiltonian_with(length, m, k)?;

    Ok(())
}
