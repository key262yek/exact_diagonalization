#[allow(unused_imports)]
use genawaiter::{sync::gen, sync_producer, yield_};
use crate::prelude::*;

pub struct PeriodicNearestXXZ{
    pub delta : f64
}

impl PeriodicNearestXXZ{
    pub fn new(delta : f64) -> Self{
        Self{
            delta,
        }
    }

    pub fn apply_to<S>(&self, state : S) -> impl Iterator<Item = (usize, f64)> + '_
        where S : State + 'static{
        gen!({
            let num = state.rep();
            let mut sum = 0f64;
            for ((i, si), (j, sj)) in state.periodic_pair_enumerate(){
                if si == sj{
                    sum -= self.delta / 2f64;
                } else {
                    sum += self.delta / 2f64;
                    let flipped = bit_flip_unsafe(num, i, j);
                    yield_!((flipped, -1f64));
                }
            }

            yield_!((num, sum));
        }).into_iter()
    }
}

pub struct PeriodicNextNearestXXZ{
    pub delta : f64
}

impl PeriodicNextNearestXXZ{
    pub fn new(delta : f64) -> Self{
        Self{
            delta,
        }
    }

    pub fn delta(&self) -> f64{
        self.delta
    }

    pub fn apply_to<S>(&self, state : S) -> impl Iterator<Item = (usize, f64)> + '_
        where S : State + 'static{
        gen!({
            let num = state.rep();
            let mut sum = 0f64;
            for ((i, si), (j, sj)) in state.periodic_pair_enumerate(){
                if si == sj{
                    sum -= self.delta / 2f64;
                } else {
                    sum += self.delta / 2f64;
                    let flipped = bit_flip_unsafe(num, i, j);
                    yield_!((flipped, -1f64));
                }
            }

            for ((i, si), (j, sj)) in state.periodic_distanced_pair_enumerate(2){
                if si != sj{
                    let flipped = bit_flip_unsafe(num, i, j);
                    yield_!((flipped, -1f64));
                }
            }

            yield_!((num, sum));
        }).into_iter()
    }
}




#[cfg(test)]
mod test {
    use super::*;
    use crate::states::SimpleState;

    #[test]
    fn test_apply_hamiltonian(){

        let delta = 2f64;
        let x: PeriodicNearestXXZ  = PeriodicNearestXXZ::new(delta);
        for data in x.apply_to(SimpleState::new(2, 4)){
            assert!((data == (1, -1f64)) || (data == (4, -1f64)) || (data == (2, 0f64)));
        }
    }
}



