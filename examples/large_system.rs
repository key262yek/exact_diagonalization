use ndarray_linalg::{Eigh, UPLO};
use exact_diagonalization::prelude::*;
use std::env;
use std::io::prelude::*;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;



fn hamiltonian_with(l : usize, m : usize, k : usize) -> Result<Array1<f64>, ()>{
    let basis_gen = BasisNK::new(EigenNumMomentum::new(m, k), l);
    let (basis, indices) = basis_gen.build().unwrap();

    let delta = 2f64;
    let xxz = PeriodicNearestXXZ::new(1f64, delta / 2f64);
    let n = basis.len();
    if n == 0{
        return Err(());
    }
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

    return hamiltonian.eigh(UPLO::Lower).map(|e| e.0).map_err(|_e| ());
}

#[allow(dead_code)]
fn hamiltonian_with_next(l : usize, m : usize, k : usize) -> Result<Array1<f64>, ()>{
    let basis_gen = BasisNK::new(EigenNumMomentum::new(m, k), l);
    let (basis, indices) = basis_gen.build().unwrap();

    let delta = 2f64;
    let xxz = PeriodicNextNearestXXZ::new(1f64, 1f64, delta / 2f64);
    let n = basis.len();
    if n == 0{
        return Err(());
    }
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

    return hamiltonian.eigh(UPLO::Lower).map(|e| e.0).map_err(|_e| ());
}



fn main() -> (){
    let args: Vec<String> = env::args().collect();

    println!("{:?}", &args);
    let length = args[1].parse::<usize>().unwrap();
    let filepath = Path::new(&args[2]);

    let output = File::create(filepath).unwrap();
    let mut writer = BufWriter::new(&output);

    for i in 0..length + 1{
        for k in 0..length{
            match hamiltonian_with(length, i, k){
                Ok(eig_v) => {
                    write!(&mut writer,"{} {} {} {:?}\n", length, i, k, eig_v).unwrap();
                },
                Err(_) => {
                    continue;
                }
            }
        }
    }
}
