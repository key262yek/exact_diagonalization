use num::integer::binomial;
use exact_diagonalization::prelude::*;


#[test]
fn test_number_conservation() -> Result<(), Error> {
    let total_ptl = 4;
    let max_state = 1 << total_ptl;
    let mut magnet_sets : Vec<Vec<SimpleState>> = Vec::new();
    let mut indices : HashMap<usize, (usize, usize)> = HashMap::new();

    let delta = 2f64;
    let t = 2f64 * delta;

    for i in 0..total_ptl + 1{
        magnet_sets.push(Vec::with_capacity(binomial(total_ptl, i)));
    }

    for n in 0..max_state {
        let state = SimpleState::new(n, total_ptl);
        let m = state.total_number();

        let idx = magnet_sets[m].len();
        magnet_sets[m].push(state);
        indices.insert(n, (m, idx));
    }

    for (m, mset) in magnet_sets.iter().enumerate(){
        let n = binomial(total_ptl, m);
        let mut hamiltonian : Array2<f64> = Array2::zeros((n, n));

        for (idx, &state) in mset.iter().enumerate(){
            for ((i, si), (j, sj)) in state.periodic_pair_enumerate() {
                if si == sj {
                    hamiltonian[[idx, idx]] -= delta / 2f64;
                } else {
                    hamiltonian[[idx, idx]] += delta / 2f64;

                    let rep2 = bit_flip(state.rep, state.length, i, j)?;
                    let idx2 = indices.get(&rep2).unwrap();
                    hamiltonian[[idx2.1, idx]] -= 1f64;
                }
            }
        }

        if m == 0 || m == 4 {
            assert_eq!(hamiltonian, arr2(&[[-t]]));
        } else if m == 1 || m == 3 {
            assert_eq!(hamiltonian, 
                arr2(&[[    0.0,   -1.0,    0.0,   -1.0],
                       [   -1.0,    0.0,   -1.0,    0.0],
                       [    0.0,   -1.0,    0.0,   -1.0],
                       [   -1.0,    0.0,   -1.0,    0.0]]));
        } else if m == 2 {
            assert_eq!(hamiltonian, 
                arr2(&[[    0.0,   -1.0,    0.0,    0.0,   -1.0,    0.0],
                       [   -1.0,      t,   -1.0,   -1.0,    0.0,   -1.0],
                       [    0.0,   -1.0,    0.0,    0.0,   -1.0,    0.0],
                       [    0.0,   -1.0,    0.0,    0.0,   -1.0,    0.0],
                       [   -1.0,    0.0,   -1.0,   -1.0,      t,   -1.0],
                       [    0.0,   -1.0,    0.0,    0.0,   -1.0,    0.0],]));
        }
    }

    return Ok(());
}

#[test]
fn test_number_conservation2() -> Result<(), Error> {
    let total_ptl = 4;
    let basis_gen = Basis::new(total_ptl);
    let (magnet_sets, indices) = basis_gen.build_n();

    let delta = 2f64;
    let t = 2f64 * delta;

    for (egn_n, mset) in magnet_sets.iter(){
        let m = egn_n.total_number();
        let n = binomial(total_ptl, m);
        let mut hamiltonian : Array2<f64> = Array2::zeros((n, n));

        for (idx, state) in mset.iter().enumerate(){
            for ((i, si), (j, sj)) in state.periodic_pair_enumerate() {
                if si == sj {
                    hamiltonian[[idx, idx]] -= delta / 2f64;
                } else {
                    hamiltonian[[idx, idx]] += delta / 2f64;

                    let rep2 = bit_flip(state.rep(), state.length, i, j)?;
                    let idx2 = indices.get(&(*egn_n, rep2)).unwrap();
                    hamiltonian[[*idx2, idx]] -= 1f64;
                }
            }
        }

        if m == 0 || m == 4 {
            assert_eq!(hamiltonian, arr2(&[[-t]]));
        } else if m == 1 || m == 3 {
            assert_eq!(hamiltonian,
                arr2(&[[    0.0,   -1.0,    0.0,   -1.0],
                       [   -1.0,    0.0,   -1.0,    0.0],
                       [    0.0,   -1.0,    0.0,   -1.0],
                       [   -1.0,    0.0,   -1.0,    0.0]]));
        } else if m == 2 {
            assert_eq!(hamiltonian,
                arr2(&[[    0.0,   -1.0,    0.0,    0.0,   -1.0,    0.0],
                       [   -1.0,      t,   -1.0,   -1.0,    0.0,   -1.0],
                       [    0.0,   -1.0,    0.0,    0.0,   -1.0,    0.0],
                       [    0.0,   -1.0,    0.0,    0.0,   -1.0,    0.0],
                       [   -1.0,    0.0,   -1.0,   -1.0,      t,   -1.0],
                       [    0.0,   -1.0,    0.0,    0.0,   -1.0,    0.0],]));
        }
    }

    return Ok(());
}

#[test]
fn test_number_conservation3() -> Result<(), Error> {
    let total_ptl = 4;
    let m = 2;
    let basis = BasisN::new(EigenNumber::new(m), total_ptl);

    let (magnet_set, indices) = basis.build();

    let delta = 2f64;
    let t = 2f64 * delta;

    let n = binomial(total_ptl, m);
    let mut hamiltonian : Array2<f64> = Array2::zeros((n, n));

    for (idx, state) in magnet_set.iter().enumerate(){
        for ((i, si), (j, sj)) in state.periodic_pair_enumerate() {
            if si == sj {
                hamiltonian[[idx, idx]] -= delta / 2f64;
            } else {
                hamiltonian[[idx, idx]] += delta / 2f64;

                let rep2 = bit_flip(state.rep(), total_ptl, i, j)?;
                if let Some(&idx2) = indices.get(&rep2){
                    hamiltonian[[idx2, idx]] -= 1f64;
                }
            }
        }
    }

    assert_eq!(hamiltonian,
        arr2(&[[    0.0,   -1.0,    0.0,    0.0,   -1.0,    0.0],
               [   -1.0,      t,   -1.0,   -1.0,    0.0,   -1.0],
               [    0.0,   -1.0,    0.0,    0.0,   -1.0,    0.0],
               [    0.0,   -1.0,    0.0,    0.0,   -1.0,    0.0],
               [   -1.0,    0.0,   -1.0,   -1.0,      t,   -1.0],
               [    0.0,   -1.0,    0.0,    0.0,   -1.0,    0.0],]));

    return Ok(());
}
