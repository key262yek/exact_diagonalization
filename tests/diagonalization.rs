use ndarray_linalg::{Eigh, UPLO, assert::close_l2};
use exact_diagonalization::prelude::*;

fn nearest_hamiltonian_with(l : usize, m : usize, k : usize) -> Result<(Array1<f64>, Array2<Complex64>), Error>{
    let basis_gen = BasisNK::new(EigenNumMomentum::new(m, k), l);
    let (basis, indices) = basis_gen.build().unwrap();

    let delta = 2f64;
    let xxz = PeriodicNearestXXZ::new(1f64, delta);
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


#[test]
fn test_diagonalization() -> Result<(), Error>{
    let t = Complex64::from(2f64.sqrt());
    let s = Complex64::from(3f64.sqrt());
    let cn = 2f64 * (6f64 - 2f64 * s).sqrt();
    let cp = 2f64 * (6f64 + 2f64 * s).sqrt();

    let length = 4;
    let m = 2;


    let (test_values, test_vectors) = nearest_hamiltonian_with(length, m, 0)?;
    let truth_values = arr1(&[2f64 - 12f64.sqrt(), 2f64 + 12f64.sqrt()]);
    let truth_vectors = arr2(&[[((-2f64 * t) / cn), ((-2f64 * t) / cp)],
                                                            [(2f64 * (1f64 - s) / cn), (2f64 * (1f64 + s) / cp)]]);
    close_l2(&test_values, &truth_values, 1e-10);
    close_l2(&test_vectors, &truth_vectors, 1e-10);
    Ok(())
}
