use exact_diagonalization::hamiltonian::{count_degeneracy_from, degeneracy_pair, prepare_energy_map};
use exact_diagonalization::states::bit_fns::period_unsafe;
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
use std::sync::{Arc, Mutex};
use std::thread;



fn rng_seed(seed : u128) -> Pcg64{
    // PCG family algorithm을 기반해서 random number를 만드는 generator를 만들어주는 함수.
    // seed : random number를 결정지을 seed.

    const INC: u128 = 0xa02bdbf7bb3c0a7ac28fa16a64abf96;
    rand_pcg::Pcg64::new(seed, INC)
}


fn normalize_factor(state : &(usize, usize)) -> Complex64{
    let p = period_unsafe(state.0, state.1);
    Complex64::from((p as f64).sqrt() / (state.1 as f64))
}

fn hamiltonian_nearest(basis : &Vec<(usize, usize)>, indices : &FnvHashMap<(EigenNumMomentum, usize), (usize, usize)>, egn_nk : EigenNumMomentum, delta : f64)
                    -> Result<Array2<Complex64>, ()>{

    let length = basis[0].1;
    let xxz = PeriodicNearestXXZ::new(1f64,  delta);
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
    let mut args = env::args().skip(1);
    // println!("{:?}", env::current_dir());
    // println!("{:?}", args);

    let interaction_info = args.next().unwrap().parse::<usize>().unwrap();
    let l = args.next().unwrap().parse::<usize>().unwrap();
    let m = args.next().unwrap().parse::<usize>().unwrap();
    let k = args.next().unwrap().parse::<usize>().unwrap();
    let delta  = args.next().unwrap().parse::<f64>().unwrap();
    let lambda = args.next().unwrap().parse::<f64>().unwrap();
    let num_en     = args.next().unwrap().parse::<usize>().unwrap();
    let org_seed = args.next().unwrap().parse::<usize>().unwrap();
    let num_thread = args.next().unwrap().parse::<usize>().unwrap();
    let en_p_thr = num_en / num_thread;
    let threshold = args.next().unwrap().parse::<f64>().unwrap();

    let p : f64 = 1f64 / (num_en as f64);
    let egn_nk = EigenNumMomentum::new(m, k);

    let total_start = Instant::now();

    let basis_gen = Basis::new(l);
    let (basis_map, indices) = basis_gen.build_light_nk();

    let filepath = format_args!("examples/output/parallel_ed_p_Index_{}_SysSize_{}_M_{}_K_{}_Delta_{:0e}_Lambda_{:0e}_Ensem_{}_Seed_{}_Thr_{:0e}.dat", interaction_info, l, m, k, delta, lambda, num_en, org_seed, threshold).to_string();
    let output = File::create(filepath).unwrap();
    let mut writer = BufWriter::new(&output);
    let mut hamiltonian_set : FnvHashMap<EigenNumMomentum, (Array2<Complex64>, Array2<Complex64>, Array1<Complex64>, Array2<Complex64>)> = FnvHashMap::default();

    let arc_basis = Arc::new(basis_map.get(&egn_nk).unwrap().clone());
    let arc_target_h0 = match interaction_info{
        0 => Arc::new(hamiltonian_nearest(&arc_basis, &indices, egn_nk, delta).unwrap()),
        1 => Arc::new(hamiltonian_next_nearest(&arc_basis, &indices, egn_nk, delta).unwrap()),
        _ => panic!(),
    };
    let (eval0, evec0) = &arc_target_h0.eigh(UPLO::Lower).unwrap();
    let min_energy_gap = eval0.iter().zip(eval0.iter().skip(1))
                            .map(|(x1, x2)| x2 - x1)
                            .min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let rtol = threshold * min_energy_gap;
    // println!("{}", min_energy_gap);

    let arc_target_h1 = Arc::new(diag_ising(&arc_basis, lambda).unwrap());
    let arc_target_cvec0 : Arc<Array2<Complex64>> = Arc::new(conjugate(evec0));
    let arc_target_ediag0 = Arc::new(eval0.map(|&x| Complex64::new(x, 0.0)));

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
        let (eval, evec) = &other_h.eigh(UPLO::Lower).unwrap();
        if count_degeneracy_from(&mut energy_map, egn, eval, rtol){
            let other_h1 = diag_ising(other_basis, lambda).unwrap();
            let cvec = conjugate(evec);
            let ediag = eval.map(|&x| Complex64::new(x, 0.0));
            hamiltonian_set.insert(egn, (other_h, other_h1, ediag, cvec));
        }
    }

    let (num_deg, map_deg) = degeneracy_pair(arc_basis.len(), &energy_map).unwrap();

    let seed = (org_seed + l + m + k + interaction_info + (delta + lambda) as usize) as u128;
    let arc_rng = Arc::new(Mutex::new(rng_seed(seed)));
    let uni = Uniform::new(-1f64, 1f64);

    let arc_result : Arc<Mutex<Array1<f64>>> = Arc::new(Mutex::new(Array1::zeros(arc_basis.len())));
    let mut handles = vec![];

    let arc_hamiltonian_set = Arc::new(hamiltonian_set);
    let arc_num_deg = Arc::new(num_deg);
    let arc_map_deg = Arc::new(map_deg);

    for thrd in 0..num_thread {
        let mutex_result = Arc::clone(&arc_result);
        let mutex_rng = Arc::clone(&arc_rng);
        let basis = Arc::clone(&arc_basis);
        let target_h0= Arc::clone(&arc_target_h0);
        let hamiltonian_set = Arc::clone(&arc_hamiltonian_set);
        let target_h1 = Arc::clone(&arc_target_h1);
        let target_cvec0 = Arc::clone(&arc_target_cvec0);
        let target_ediag0 = Arc::clone(&arc_target_ediag0);
        let num_deg = Arc::clone(&arc_num_deg);
        let map_deg = Arc::clone(&arc_map_deg);
        let handle = thread::spawn( move ||{
            let start = Instant::now();
            for _en in 0..en_p_thr{
                let r = {
                    let mut rng = mutex_rng.lock().unwrap();
                    rng.sample(uni)
                };
                let mut total : Array1<Complex64> = Array1::zeros(basis.len());

                let change = compute_energy_change(&target_h0, &target_h1, r, &target_ediag0, &target_cvec0);
                total = total + &change * &(*num_deg);

                match map_deg.get(&egn_nk){
                    Some(v) => {
                        for (from, to) in v{
                            total[*to]+= change[*from] * num_deg[*to];
                        }
                    },
                    None => {},
                }

                for (egn, (h0, h1, ediag, cvec)) in hamiltonian_set.iter(){
                    let change = compute_energy_change(h0, h1, r, ediag, cvec);
                    for (from, to) in map_deg.get(egn).unwrap(){
                        total[*to]+= change[*from] * num_deg[*to];
                    }
                }

                {
                    let mut result = mutex_result.lock().unwrap();
                    for (res, &work) in result.iter_mut().zip(total.map(|x| if x.re > rtol {p} else {0.0}).iter()){
                        *res += work;
                    }
                }
            }
            println!("Thread {} works with ensemble {} in time {:?}", thrd, en_p_thr, start.elapsed());
        });
        handles.push(handle);
    };
    for handle in handles{
        handle.join().unwrap();
    }

    let result = arc_result.lock().unwrap();
    for i in 0..arc_basis.len(){
        let v = eval0[i];
        let pv = result[i];
        let d = arc_num_deg[i];
        if d < 0f64{
            continue;
        }
        write!(&mut writer,"{:.05e}\t{:.05e}\n", v, pv).unwrap();
    }
    println!("Total time : {:?}", total_start.elapsed());
}

// fn main(){
//     println!("Hello World");
// }
