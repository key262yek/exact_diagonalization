use exact_diagonalization::states::bit_fns::bit_flip_unsafe;
use criterion::{criterion_group, criterion_main, Criterion};
use exact_diagonalization::{prelude::*};
use ndarray_linalg::eigh::EighInplace;


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

fn use_eigh_inplace(h : &mut Array2<Complex64>){
    // Collecting 100 samples in estimated 35.895 s (100 iterations)
    //  time:   [350.43 ms 352.22 ms 354.30 ms]
    h.eigh_inplace(UPLO::Lower).unwrap();
}

fn bench_eigh_inplace(c: &mut Criterion) {
    let mut hamiltonian = hamiltonian_with(15, 7, 0).unwrap();
    c.bench_function("use eigh_inplace", |b| {
        b.iter(|| use_eigh_inplace(&mut hamiltonian))
    });
}

fn use_eigh(h : &Array2<Complex64>){
    // Collecting 100 samples in estimated 28.429 s (100 iterations)
    // time:   [285.95 ms 288.62 ms 292.02 ms]
    //
    h.eigh(UPLO::Lower).unwrap();
}

fn bench_eigh(c: &mut Criterion) {
    let hamiltonian = hamiltonian_with(15, 7, 0).unwrap();
    c.bench_function("use eigh", |b| {
        b.iter(|| use_eigh(&hamiltonian))
    });
}

criterion_group!(benches, bench_eigh, bench_eigh_inplace);
criterion_main!(benches);
