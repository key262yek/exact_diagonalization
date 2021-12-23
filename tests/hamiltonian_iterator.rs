use exact_diagonalization::states::bit_fns::bit_flip_unsafe;
use genawaiter::{sync::gen, yield_};
use exact_diagonalization::prelude::*;


#[test]
fn test_genawaiter() -> Result<(), Error> {
    let total_ptl = 4;
    let max_state = 1 << total_ptl;
    let mut hamiltonian = Array2::<f64>::zeros((max_state, max_state));
    let delta = 2f64;
    let t = 2f64 * delta;

    for state in 0..max_state {
        let iterator = SimpleState::new(state, total_ptl);

        let hamilton_iterator = gen!({
            let mut sum = 0f64;
            for ((i, si), (j, sj)) in iterator.periodic_pair_enumerate(){
                if si == sj{
                    sum -= delta / 2f64;
                } else {
                    sum += delta / 2f64;
                    let flipped = bit_flip_unsafe(state, i, j);
                    yield_!((flipped, -1f64));
                }
            }

            yield_!((state, sum));
        });

        for (state2, val) in hamilton_iterator{
            hamiltonian[[state2, state]] += val;
        }
    }

    assert_eq!(hamiltonian,
        arr2(&[[   -t,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
               [  0.0, -1.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0, -1.0,  0.0,    t, -1.0,  0.0,  0.0, -1.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0, -1.0,  0.0],
               [  0.0, -1.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0, -1.0,  0.0,  0.0, -1.0,  0.0,  0.0, -1.0,    t,  0.0, -1.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0, -1.0,  0.0],
               [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   -t]]));
    return Ok(());
}

#[test]
fn test_genawaiter2() -> Result<(), Error> {
    let total_ptl = 4;
    let max_state = 1 << total_ptl;
    let mut hamiltonian = Array2::<f64>::zeros((max_state, max_state));
    let delta = 2f64;
    let t = 2f64 * delta;
    let xxz = PeriodicNearestXXZ::new(delta);

    for idx in 0..max_state {
        let state = SimpleState::new(idx, total_ptl);

        for (idx2, val) in xxz.apply_to(&state){
            hamiltonian[[idx2, idx]] += val;
        }
    }

    assert_eq!(hamiltonian,
        arr2(&[[   -t,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
               [  0.0, -1.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0, -1.0,  0.0,    t, -1.0,  0.0,  0.0, -1.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0, -1.0,  0.0],
               [  0.0, -1.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0, -1.0,  0.0,  0.0, -1.0,  0.0,  0.0, -1.0,    t,  0.0, -1.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0, -1.0,  0.0],
               [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0],
               [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   -t]]));
    return Ok(());
}