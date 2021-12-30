use exact_diagonalization::hamiltonian::{count_degeneracy_from, degeneracy_pair, prepare_energy_map};
use exact_diagonalization::states::bit_fns::period_unsafe;
use fnv::FnvHashMap;
use exact_diagonalization::prelude::*;
use std::env;
use std::io::prelude::*;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;



fn normalize_factor(state : &(usize, usize)) -> Complex64{
    let p = period_unsafe(state.0, state.1);
    Complex64::from((p as f64).sqrt() / (state.1 as f64))
}

fn hamiltonian_with(basis : &Vec<(usize, usize)>, indices : &FnvHashMap<(EigenNumMomentum, usize), (usize, usize)>, egn_nk : EigenNumMomentum, delta : f64)
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

fn diag_ising(basis : &Vec<(usize, usize)>, delta : f64) -> Result<Array2<Complex64>, ()>{
    let ising = PeriodicIsing::new(delta);
    let n = basis.len();
    let mut hamiltonian = Array2::zeros((n, n));
    for (idx, state) in basis.iter().enumerate(){
        let value = ising.apply_to(state);
        hamiltonian[[idx, idx]] += Complex64::from(value);
    }

    return Ok(hamiltonian);
}

fn compute_energy_change(h0 : &Array2<Complex64>, h1 : &Array2<Complex64>, r : f64, ediag : &Array1<Complex64>, cvec : &Array2<Complex64>) -> Array1<Complex64>{

    let h2 = (h1 * r) + h0;
    let (eval2, evec2) = &h2.eigh(UPLO::Lower).unwrap();
    let unitary1 = Array2::from_diag(&eval2.map(|&x| Complex64::new(0.0, x).exp()));
    let x1 = cvec.dot(evec2);
    let x2 : Array2<Complex64> = conjugate(&x1);
    let right = x1.dot(&unitary1).dot(&x2);
    let right_abs = right.map(|x| x * x.conj());

    right_abs.apply(ediag) - ediag
}

fn main() -> (){
    let args: Vec<String> = env::args().collect();
    // println!("{:?}", env::current_dir());
    // println!("{:?}", args);

    let l = args[1].parse::<usize>().unwrap();
    let m = args[2].parse::<usize>().unwrap();
    let k = args[3].parse::<usize>().unwrap();
    let delta  = args[4].parse::<f64>().unwrap();
    let lambda = args[5].parse::<f64>().unwrap();
    let r = args[6].parse::<f64>().unwrap();
    let filepath = Path::new(&args[7]);

    let rtol : f64 = 1e-4;
    let p : f64 = 1f64;
    let egn_nk = EigenNumMomentum::new(m, k);


    let basis_gen = Basis::new(l);
    let (basis_map, indices) = basis_gen.build_light_nk();

    let output = File::create(filepath).unwrap();
    let mut writer = BufWriter::new(&output);
    let mut hamiltonian_set : FnvHashMap<EigenNumMomentum, (Array2<Complex64>, Array2<Complex64>, Array1<Complex64>, Array2<Complex64>)> = FnvHashMap::default();

    let basis = basis_map.get(&egn_nk).unwrap();
    let target_h0 = hamiltonian_with(basis, &indices, egn_nk, delta).unwrap();
    let (eval0, evec0) = &target_h0.eigh(UPLO::Lower).unwrap();
    let target_h1 = diag_ising(basis, lambda).unwrap();
    let target_cvec0 : Array2<Complex64> = conjugate(evec0);
    let target_ediag0 = eval0.map(|&x| Complex64::new(x, 0.0));

    let mut energy_map = prepare_energy_map(egn_nk,eval0, rtol);

    for (&egn, other_basis) in basis_map.iter(){
        if egn == egn_nk{
            continue;
        }

        let other_h = hamiltonian_with(other_basis, &indices, egn, delta).unwrap();
        let (eval, evec) = &other_h.eigh(UPLO::Lower).unwrap();
        if count_degeneracy_from(&mut energy_map, egn, eval, rtol){
            let other_h1 = diag_ising(other_basis, lambda).unwrap();
            let cvec = conjugate(evec);
            let ediag = eval.map(|&x| Complex64::new(x, 0.0));
            hamiltonian_set.insert(egn, (other_h, other_h1, ediag, cvec));
        }
    }

    let (num_deg, map_deg) = degeneracy_pair(basis.len(), &energy_map).unwrap();

    let mut total : Array1<Complex64> = Array1::zeros(basis.len());

    let change0 = compute_energy_change(&target_h0, &target_h1, r, &target_ediag0, &target_cvec0);
    println!("{:?}\n{:?}\n\n", egn_nk, change0.map(|x| x.re));
    total = total + &change0 * &num_deg;

    match map_deg.get(&egn_nk){
        Some(v) => {
            for (from, to) in v{
                total[*to]+= change0[*from] * num_deg[*to];
            }
        },
        None => {},
    }

    for (egn, (h0, h1, ediag, cvec)) in hamiltonian_set.iter(){
        let change = compute_energy_change(h0, h1, r, ediag, cvec);
        for (from, to) in map_deg.get(egn).unwrap(){
            total[*to]+= change[*from] * num_deg[*to];
        }
        println!("{:?}\n{:?}\n\n", egn, change.map(|x| x.re));
    }

    for i in 0..basis.len(){
        let e = eval0[i];
        let work = total[i].re;
        let c = change0[i].re;
        let d = num_deg[i];
        if d < 0.0f64{
            continue;
        }
        write!(&mut writer,"{:.05e}\t{:.05e}\t{:.5e}\n", e, work, c).unwrap();
    }
    // for (energy, work) in eval0.iter()
    //                                 .zip(total.iter()
    //                                   .zip(num_deg.iter())
    //                                   .filter(|(_x, d)| **d > 0.0)
    //                                   .map(|(x, _d)| x.re)){
    //     write!(&mut writer,"{:.05e}\t{:.05e}\n", energy, work).unwrap();
    // }
}
