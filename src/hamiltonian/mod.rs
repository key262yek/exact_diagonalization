#[allow(unused_imports)]
use genawaiter::{sync::gen, sync_producer, yield_};
use crate::prelude::*;

#[derive(Copy, Clone, Debug)]
pub struct PeriodicIsing{
    pub delta : f64,
}

impl PeriodicIsing{
    pub fn new(delta : f64) -> Self{
        Self{
            delta,
        }
    }

    pub fn apply_to<S>(&self, state : &S) -> f64
        where S : State{

        let mut sum = 0f64;
        for (si, sj) in state.periodic_pair_iter(){
            if si == sj{
                sum -= self.delta / 2f64;
            } else {
                sum += self.delta / 2f64;
            }
        }

        return sum;
    }
}

#[derive(Copy, Clone, Debug)]
pub struct PeriodicNearestXXZ{
    pub delta_x : f64,
    pub delta_z : f64,
}

impl PeriodicNearestXXZ{
    pub fn new(delta_x : f64, delta_z : f64) -> Self{
        Self{
            delta_x,
            delta_z,
        }
    }

    pub fn apply_to<'a, S>(&'a self, state : &'a S) -> impl Iterator<Item = (usize, f64)> + 'a
        where S : State + 'static{
        gen!({
            let num = state.rep();
            let mut sum = 0f64;
            for ((i, si), (j, sj)) in state.periodic_pair_enumerate(){
                if si == sj{
                    sum -= self.delta_z / 2f64;
                } else {
                    sum += self.delta_z / 2f64;
                    let flipped = bit_flip_unsafe(num, i, j);
                    yield_!((flipped, -self.delta_x));
                }
            }

            yield_!((num, sum));
        }).into_iter()
    }
}

#[derive(Copy, Clone, Debug)]
pub struct PeriodicNextNearestXXZ{
    pub delta_x1 : f64,
    pub delta_x2 : f64,
    pub delta_z : f64,
}

impl PeriodicNextNearestXXZ{
    pub fn new(delta_x1 : f64, delta_x2 : f64, delta_z : f64) -> Self{
        Self{
            delta_x1,
            delta_x2,
            delta_z,
        }
    }

    pub fn apply_to<'a, S>(&'a self, state : &'a S) -> impl Iterator<Item = (usize, f64)> + 'a
        where S : State + 'static{
        gen!({
            let num = state.rep();
            let mut sum = 0f64;
            for ((i, si), (j, sj)) in state.periodic_pair_enumerate(){
                if si == sj{
                    sum -= self.delta_z / 2f64;
                } else {
                    sum += self.delta_z / 2f64;
                    let flipped = bit_flip_unsafe(num, i, j);
                    yield_!((flipped, -self.delta_x1));
                }
            }

            for ((i, si), (j, sj)) in state.periodic_distanced_pair_enumerate(2){
                if si != sj{
                    let flipped = bit_flip_unsafe(num, i, j);
                    yield_!((flipped, -self.delta_x2));
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
        let x: PeriodicNearestXXZ  = PeriodicNearestXXZ::new(1f64, delta);
        for data in x.apply_to(&SimpleState::new(2, 4)){
            assert!((data == (1, -1f64)) || (data == (4, -1f64)) || (data == (2, 0f64)));
        }
    }
}



