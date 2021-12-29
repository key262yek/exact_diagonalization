use fnv::FnvHashMap;
use ndarray_linalg::close_l2;
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
fn test_quench() -> (){
    let l = 4;
    let m = 2;
    let k = 0;
    let delta  = 2f64;
    let lambda = 1f64;
    let rtol = 1e-6;

    let basis_gen = BasisNK::new(EigenNumMomentum::new(m, k), l);
    let (basis, indices) = basis_gen.build().unwrap();

    let h0 = hamiltonian_with(&basis, &indices, delta).unwrap();
    let (eval0, evec0) = &h0.eigh(UPLO::Lower).unwrap();
    let conj_evec0 : Array2<Complex64> = conjugate(evec0);
    let energy_diag = Array2::from_diag(&eval0.map(|&x| Complex64::new(x, 0.0)));

    let h1 = diag_ising(&basis, lambda).unwrap();

    let mut rng = rng_seed(124178297891748914u128);
    let uni = Uniform::new(-1f64, 1f64);
    let r = rng.sample(uni);
    let h2 = (&h1 * r) + &h0;
    let truth_h2 : Array2<Complex64> = arr2(&[
        [Complex64 { re: 0.0, im: 0.0 }, Complex64 { re: -2.8284271247461903, im: 0.0 }],
        [Complex64 { re: -2.82842712474619, im: 0.0 }, Complex64 { re: 3.6085624472639157, im: 0.0 }]
    ]);
    close_l2(&h2, &truth_h2, rtol);

    let (eval2, evec2) = &h2.eigh(UPLO::Lower).unwrap();
    let truth_eval2 : Array1<f64> = arr1(&[-1.55063021095925, 5.159192658223166]);
    let truth_evec2 : Array2<Complex64> = arr2(&[
        [Complex64 { re: -0.8768702689733415, im: 0.0 }, Complex64 { re: -0.48072708618364707, im: 0.0 }],
        [Complex64 { re: -0.48072708618364707, im: 0.0 }, Complex64 { re: 0.8768702689733415, im: 0.0 }]
    ]);
    close_l2(eval2, &truth_eval2, rtol);
    close_l2(evec2, &truth_evec2, rtol);

    let unitary1 = Array2::from_diag(&eval2.map(|&x| Complex64::new(0.0, x).exp()));
    let truth_uni : Array2<Complex64> = arr2(&[
        [Complex64 { re: 0.020164749030229835, im: -0.9997966707768874 }, Complex64 { re: 0.0, im: 0.0 }],
        [Complex64 { re: 0.0, im: 0.0 }, Complex64 { re: 0.4320851980779638, im: -0.9018327902676453 }]
    ]);
    close_l2(&unitary1, &truth_uni, rtol);

    let x1 = conj_evec0.dot(evec2);
    let truth_x1 : Array2<Complex64> = arr2(&[
        [Complex64 { re: 0.9997161886224049, im: 0.0 }, Complex64 { re: 0.023823144341003966, im: 0.0 }],
        [Complex64 { re: -0.023823144341003966, im: 0.0 }, Complex64 { re: 0.9997161886224049, im: 0.0 }]
    ]);
    close_l2(&x1, &truth_x1, rtol);

    let x2 : Array2<Complex64> = conjugate(&x1);
    let truth_x2 : Array2<Complex64> = arr2(&[
        [Complex64 { re: 0.9997161886224049, im: -0.0 }, Complex64 { re: -0.023823144341003966, im: -0.0 }],
        [Complex64 { re: 0.023823144341003966, im: -0.0 }, Complex64 { re: 0.9997161886224049, im: -0.0 }]
    ]);
    close_l2(&x2, &truth_x2, rtol);

    let right = x1.dot(&unitary1).dot(&x2);
    let truth_right : Array2<Complex64> = arr2(&[
        [Complex64 { re: 0.020398531270699308, im: -0.9997410721400063 }, Complex64 { re: 0.009810455205422956, im: 0.002333145304407931 }],
        [Complex64 { re: 0.009810455205422956, im: 0.002333145304407931 }, Complex64 { re: 0.43185141583749437, im: -0.9018883889045265 }]
    ]);
    close_l2(&right, &truth_right, rtol);

    let left : Array2<Complex64> = conjugate(&right);
    let truth_left : Array2<Complex64> = arr2(&[
        [Complex64 { re: 0.020398531270699308, im: 0.9997410721400063 }, Complex64 { re: 0.009810455205422956, im: -0.002333145304407931 }],
        [Complex64 { re: 0.009810455205422956, im: -0.002333145304407931 }, Complex64 { re: 0.43185141583749437, im: 0.9018883889045265 }]
    ]);
    close_l2(&left, &truth_left, rtol);

    let change = (left.dot(&energy_diag).dot(&right) - &energy_diag).into_diag();
    let truth_change : Array1<Complex64> = arr1(&[
        Complex64 { re: 0.0007045192755643637, im: 0.000000000000000003469446951953614 },
        Complex64 { re: -0.0007045192755645857, im: 0.0 }]);
    close_l2(&change, &truth_change, rtol);

}
