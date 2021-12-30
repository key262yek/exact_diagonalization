use ndarray_linalg::Norm;
use ndarray_linalg::LinearOperator;
use fnv::FnvHashMap;


use ndarray_linalg::{Eigh, UPLO, generate::conjugate};
use exact_diagonalization::prelude::*;
use rand_pcg::Pcg64;
use rand::distributions::Uniform;
use rand::Rng;




fn rng_seed(seed : u128) -> Pcg64{
    // PCG family algorithm을 기반해서 random number를 만드는 generator를 만들어주는 함수.
    // seed : random number를 결정지을 seed.

    const INC: u128 = 0xa02bdbf7bb3c0a7ac28fa16a64abf96;
    rand_pcg::Pcg64::new(seed, INC)
}


fn hamiltonian_with(basis : &Vec<EigenState<EigenNumMomentum>>, indices : &FnvHashMap<usize, (usize, usize)>, delta : f64)
                    -> Result<Array2<Complex64>, ()>{

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

    return Ok(hamiltonian);
}

fn diag_ising(basis : &Vec<EigenState<EigenNumMomentum>>, delta : f64)
                    -> Result<Array2<Complex64>, ()>{
    let ising = PeriodicIsing::new(delta);
    let n = basis.len();
    let mut hamiltonian = Array2::zeros((n, n));
    for (idx, state) in basis.iter().enumerate(){
        let value = ising.apply_to(state);
        hamiltonian[[idx, idx]] += Complex64::from(value);
    }

    return Ok(hamiltonian);
}

#[test]
fn test_quench_faster() -> (){
    let l = 10;
    let m = 4;
    let k = 0;
    let delta  = 2f64;
    let lambda = 1f64;
    let num_en = 300;
    let rtol = 1e-2;
    let p : f64 = 1f64 / (num_en as f64);

    let basis_gen = BasisNK::new(EigenNumMomentum::new(m, k), l);
    let (basis, indices) = basis_gen.build().unwrap();

    let h0 = hamiltonian_with(&basis, &indices, delta).unwrap();
    let (eval0, evec0) = &h0.eigh(UPLO::Lower).unwrap();
    let conj_evec0 : Array2<Complex64> = conjugate(evec0);
    let energy_diag = eval0.map(|&x| Complex64::new(x, 0.0));

    let h1 = diag_ising(&basis, lambda).unwrap();

    let mut rng = rng_seed(124178297891748914u128);
    let uni = Uniform::new(-1f64, 1f64);

    let mut result : Array1<f64> = Array1::zeros(basis.len());
    for _ in 0..num_en{

        let r = rng.sample(uni);
        let h2 = (&h1 * r) + &h0;

        let (eval2, evec2) = &h2.eigh(UPLO::Lower).unwrap();

        let unitary1 = Array2::from_diag(&eval2.map(|&x| Complex64::new(0.0, x).exp()));
        let x1 = conj_evec0.dot(evec2);
        let x2 : Array2<Complex64> = conjugate(&x1);
        let right = x1.dot(&unitary1).dot(&x2);
        let right_abs = right.map(|x| x * x.conj());

        let change = right_abs.apply(&energy_diag) - &energy_diag;

        result = result + change.map(|x| if x.re > 0.000001 {p} else {0.0});
    }

    let truth_res : Array1<f64> = arr1(&[0.9985, 0.9983, 0.9976, 0.6061, 0.9939, 0.9945,  0.0, 0.9972, 0.9972, 0.8416, 0.34, 0.0, 0.0, 0.9983, 0.0, 0.0, 0.441, 0.0, 0.0,  0.9987, 0.0, 0.0]);
    let truth_res2 : Array1<f64> = arr1(&[0.9985, 0.9983, 0.9976, 0.9939, 0.6061, 0.9945, 0.0, 0.9972, 0.9972, 0.8416, 0.0, 0.34, 0.0, 0.9983, 0.0, 0.0, 0.0, 0.441, 0.0, 0.9987, 0.0, 0.0]);

    let tol = (&result - &truth_res).norm_l2() / &truth_res.norm_l2();
    let tol2 = (&result - &truth_res2).norm_l2() / &truth_res2.norm_l2();
    assert!((tol < rtol) || (tol2 < rtol));
}
