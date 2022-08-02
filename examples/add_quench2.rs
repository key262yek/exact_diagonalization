use fnv::FnvHashMap;
use exact_diagonalization::prelude::*;
use std::env;
use std::io::prelude::*;
use std::fs::File;
use std::io::BufWriter;
use rand_pcg::Pcg64;
use rand::distributions::Uniform;
use rand::Rng;
use std::time::Instant;



fn rng_seed(seed : u128) -> Pcg64{
    // PCG family algorithm을 기반해서 random number를 만드는 generator를 만들어주는 함수.
    // seed : random number를 결정지을 seed.

    const INC: u128 = 0xa02bdbf7bb3c0a7ac28fa16a64abf96;
    rand_pcg::Pcg64::new(seed, INC)
}


fn hamiltonian_nearest(basis : &Vec<NumMomentumState>, indices : &FnvHashMap<Representation<EigenNumMomentum>, (usize, usize)>, delta : f64)
                    -> Result<Array2<Complex64>, ()>{

    let xxz = PeriodicNearestXXZ::new(1f64, delta);
    let omega_k = basis[0].phase_factor();

    let n = basis.len();
    let eigen_v = basis[0].value();
    let mut hamiltonian : Array2<Complex64> = Array2::zeros((n, n));

    for (idx, state) in basis.iter().enumerate(){
        let normal_f1 = state.normalize_factor();
        for (rep2, value) in xxz.apply_to(state){
            if let Some((idx2, d)) = indices.get(&Representation(eigen_v, rep2)){
                let normal_f2 = basis[*idx2].normalize_factor();
                hamiltonian[[*idx2, idx]] += Complex64::from(value) * normal_f1 / normal_f2 * omega_k.powu(*d as u32);
            }
        }
    }

    return Ok(hamiltonian);
}


fn hamiltonian_next_nearest(basis : &Vec<NumMomentumState>, indices : &FnvHashMap<Representation<EigenNumMomentum>, (usize, usize)>, delta : f64)
                    -> Result<Array2<Complex64>, ()>{

    let xxz = PeriodicNextNearestXXZ::new(1f64, 1f64, delta);
    let omega_k = basis[0].phase_factor();

    let n = basis.len();
    let eigen_v = basis[0].value();
    let mut hamiltonian : Array2<Complex64> = Array2::zeros((n, n));

    for (idx, state) in basis.iter().enumerate(){
        let normal_f1 = state.normalize_factor();
        for (rep2, value) in xxz.apply_to(state){
            if let Some((idx2, d)) = indices.get(&Representation(eigen_v, rep2)){
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

fn main() -> (){
    let mut args = env::args().skip(1);

    if args.len() != 9{
        println!("InteractionInfo NumParticle M K Delta Lambda N_Ensemble Seed Threshold");
        panic!();
    }

    let interaction_info = args.next().unwrap().parse::<usize>().unwrap();
    let l = args.next().unwrap().parse::<usize>().unwrap();
    let m = args.next().unwrap().parse::<usize>().unwrap();
    let k = args.next().unwrap().parse::<usize>().unwrap();
    let delta  = args.next().unwrap().parse::<f64>().unwrap();
    let lambda = args.next().unwrap().parse::<f64>().unwrap();
    let num_en     = args.next().unwrap().parse::<usize>().unwrap();
    let org_seed = args.next().unwrap().parse::<usize>().unwrap();
    let threshold = args.next().unwrap().parse::<f64>().unwrap();

    let p : f64 = 1f64 / (num_en as f64);
    let egn_nk = EigenNumMomentum::new(m, k);

    let filepath = format_args!("examples/output/ed_pure_p_Index_{}_SysSize_{}_M_{}_K_{}_Delta_{:0e}_Lambda_{:0e}_Ensem_{}_Seed_{}_Thr_{:0e}.dat", interaction_info, l, m, k, delta, lambda, num_en, org_seed, threshold).to_string();
    let output = File::create(filepath).unwrap();
    let mut writer = BufWriter::new(&output);


    let basis_gen = BasisNK::new(egn_nk, l);
    let (basis, indices) = basis_gen.build().unwrap();


    let h0 = match interaction_info{
        0 => hamiltonian_nearest(&basis, &indices, delta).unwrap(),
        1 => hamiltonian_next_nearest(&basis, &indices, delta).unwrap(),
        _ => panic!(),
    };

    let (eval0, evec0) = &h0.eigh(UPLO::Lower).unwrap();
    let conj_evec0 : Array2<Complex64> = conjugate(evec0);
    let energy_diag = eval0.map(|&x| Complex64::new(x, 0.0));

    // let min_energy_gap = eval0.iter().zip(eval0.iter().skip(1))
    //                         .map(|(x1, x2)| if (x2 - x1).abs() < 1e-9 {std::f64::MAX} else {x2 - x1})
    //                         .min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    // let rtol = (threshold * min_energy_gap).max(1e-4);
    let rtol = 1e-6;

    let h1 = diag_ising(&basis, lambda).unwrap();

    let mut rng = rng_seed(124178297891748914u128);
    let uni = Uniform::new(-1f64, 1f64);

    let mut result : Array1<f64> = Array1::zeros(basis.len());
    for _ in 0..num_en{
        let start = Instant::now();

        let r = rng.sample(uni);
        let h2 = (&h1 * r) + &h0;

        let (eval2, evec2) = &h2.eigh(UPLO::Lower).unwrap();

        let unitary1 = Array2::from_diag(&eval2.map(|&x| Complex64::new(0.0, x).exp()));
        let x1 = conj_evec0.dot(evec2);
        let x2 : Array2<Complex64> = conjugate(&x1);
        let right = x1.dot(&unitary1).dot(&x2);
        let right_abs = right.map(|x| x * x.conj());

        let change = right_abs.apply(&energy_diag) - &energy_diag;

        result = result + change.map(|x| if x.re > rtol {p} else {0.0});

        println!("{:?}", start.elapsed());
    }

    for (v, pv) in eval0.iter().zip(result.iter()){
        write!(&mut writer,"{:.05e}\t{:.05e}\n", v, pv).unwrap();
    }
}
