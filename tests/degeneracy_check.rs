use ndarray_linalg::eigh::EigValsh;
use exact_diagonalization::hamiltonian::{count_degeneracy_from, prepare_energy_map};
use exact_diagonalization::states::bit_fns::period_unsafe;
use fnv::FnvHashMap;
use exact_diagonalization::prelude::*;










fn normalize_factor(state : &(usize, usize)) -> Complex64{
    let p = period_unsafe(state.0, state.1);
    Complex64::from((p as f64).sqrt() / (state.1 as f64))
}

fn hamiltonian_nearest(basis : &Vec<(usize, usize)>, indices : &FnvHashMap<(EigenNumMomentum, usize), (usize, usize)>, egn_nk : EigenNumMomentum, delta : f64)
                    -> Result<Array2<Complex64>, ()>{

    let length = basis[0].1;
    let xxz = PeriodicNearestXXZ::new(1f64, delta);
    let omega_k = egn_nk.phase_factor(length);

    let n = basis.len();
    let mut hamiltonian : Array2<Complex64> = Array2::zeros((n, n));

    for (idx, state) in basis.iter().enumerate(){
        let normal_f1 = normalize_factor(state);
        for (rep2, value) in xxz.apply_to(state){
            if let Some((idx2, d)) = indices.get(&(egn_nk, rep2)){
                let state2 = basis[*idx2];
                let normal_f2 = normalize_factor(&state2);
                hamiltonian[[*idx2, idx]] += Complex64::from(value) * normal_f1 / normal_f2 * omega_k.powu(*d as u32);
            }
        }
    }

    return Ok(hamiltonian);
}

fn hamiltonian_next_nearest(basis : &Vec<(usize, usize)>, indices : &FnvHashMap<(EigenNumMomentum, usize), (usize, usize)>, egn_nk : EigenNumMomentum, delta : f64)
                    -> Result<Array2<Complex64>, ()>{

    let length = basis[0].1;
    let xxz = PeriodicNextNearestXXZ::new(1f64, 1f64, delta);
    let omega_k = egn_nk.phase_factor(length);

    let n = basis.len();
    let mut hamiltonian : Array2<Complex64> = Array2::zeros((n, n));

    for (idx, state) in basis.iter().enumerate(){
        let normal_f1 = normalize_factor(state);
        for (rep2, value) in xxz.apply_to(state){
            if let Some((idx2, d)) = indices.get(&(egn_nk, rep2)){
                let state2 = basis[*idx2];
                let normal_f2 = normalize_factor(&state2);
                hamiltonian[[*idx2, idx]] += Complex64::from(value) * normal_f1 / normal_f2 * omega_k.powu(*d as u32);
            }
        }
    }

    return Ok(hamiltonian);
}

#[test]
#[ignore]
fn test_degeneracy_unit() -> (){
    let interaction_info = 0;
    let l = 15;
    let m = 4;
    let k = 1;
    let delta  = 2.0;
    let threshold = 1e-1;


    let egn_nk = EigenNumMomentum::new(m, k);

    let basis_gen = Basis::new(l);
    let (basis_map, indices) = basis_gen.build_light_nk();

    let basis = basis_map.get(&egn_nk).unwrap();
    let target_h0 = match interaction_info{
        0 => hamiltonian_nearest(basis, &indices, egn_nk, delta).unwrap(),
        1 => hamiltonian_next_nearest(basis, &indices, egn_nk, delta).unwrap(),
        _ => panic!(),
    };
    let eval0 = &target_h0.eigvalsh(UPLO::Lower).unwrap();
    for (i, &v) in eval0.iter().enumerate(){
        println!("{} {}", i, v);
    }
    println!("");
    let min_energy_gap = eval0.iter().zip(eval0.iter().skip(1))
                            .map(|(x1, x2)| x2 - x1)
                            .min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let rtol : f64 = (threshold * min_energy_gap).max(1e-5);

    let mut energy_map = prepare_energy_map(egn_nk,eval0, rtol);

    for (&egn, other_basis) in basis_map.iter(){
        if egn == egn_nk{
            continue;
        }

        let other_h = match interaction_info{
            0 => hamiltonian_nearest(other_basis, &indices, egn, delta).unwrap(),
            1 => hamiltonian_next_nearest(other_basis, &indices, egn, delta).unwrap(),
            _ => panic!(),
        };
        let eval= &other_h.eigvalsh(UPLO::Lower).unwrap();
        if count_degeneracy_from(&mut energy_map, egn, eval, rtol){
        }
    }

    for vector in energy_map.values(){
        println!("{:?}", vector);
    }
}
