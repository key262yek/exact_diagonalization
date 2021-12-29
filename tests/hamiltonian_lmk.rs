use num_traits::identities::Zero;
use ndarray_linalg::assert::close_l2;
use exact_diagonalization::{
    states::bit_fns::bit_flip_unsafe,
    prelude::*};

fn hamiltonian_with(l : usize, m : usize, k : usize) -> Result<Array2<Complex64>, Error>{
    let basis_gen = BasisNK::new(EigenNumMomentum::new(m, k), l);
    let (basis, indices) = basis_gen.build().unwrap();

    let delta = Complex64::from(2f64);
    let omega_k = basis[0].phase_factor();

    let n = basis.len();
    let mut hamiltonian : Array2<Complex64> = Array2::zeros((n, n));

    for (idx, state) in basis.iter().enumerate(){
        let normal_f1 = state.normalize_factor();
        for ((i, si), (j, sj)) in state.periodic_pair_enumerate() {
            if si == sj {
                hamiltonian[[idx, idx]] -= delta / 2f64;
            } else {
                hamiltonian[[idx, idx]] += delta / 2f64;

                let rep2 = bit_flip_unsafe(state.rep(), i, j);
                if let Some((idx2, d)) = indices.get(&rep2){
                    let normal_f2 = basis[*idx2].normalize_factor();
                    hamiltonian[[*idx2, idx]] += -normal_f1 / normal_f2 * omega_k.powu(*d as u32);
                }
            }
        }
    }

    return Ok(hamiltonian)
}

#[test]
fn test_hamiltonian_lmk() -> Result<(), Error>{
    let delta = Complex64::from(2f64);
    let t = Complex64::from(2f64 * delta);
    let s = Complex64::from(8f64.sqrt());
    let z = Complex64::zero();

    let length = 4;
    let m = 2;

    let test_h0 = hamiltonian_with(length, m, 0)?;
    let truth_h0 = arr2(&[[z, -s],
                                                        [-s, t]]);
    close_l2(&test_h0, &truth_h0, 1e-10);

    let test_h1 = hamiltonian_with(length, m, 1)?;
    let truth_h1 = arr2(&[[z]]);
    close_l2(&test_h1, &truth_h1, 1e-10);

    let test_h2 = hamiltonian_with(length, m, 2)?;
    let truth_h2 = arr2(&[[z, z],
                                                        [z, t]]);
    close_l2(&test_h2, &truth_h2, 1e-10);

    let test_h3 = hamiltonian_with(length, m, 3)?;
    let truth_h3 = arr2(&[[z]]);
    close_l2(&test_h3, &truth_h3, 1e-10);
    Ok(())
}
