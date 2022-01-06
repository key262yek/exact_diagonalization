#[allow(unused_imports)]
use std::hash::Hash;
use genawaiter::{sync::gen, yield_};
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

    pub fn apply_to<S, T>(&self, state : &S) -> f64
        where S : State<T>,
              T : EigenValue{

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

    pub fn apply_to<'a, S, T>(&'a self, state : &'a S) -> impl Iterator<Item = (usize, f64)> + 'a
        where S : State<T>,
              T : EigenValue{
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

    pub fn apply_to<'a, S, T>(&'a self, state : &'a S) -> impl Iterator<Item = (usize, f64)> + 'a
        where S : State<T>,
              T : EigenValue{
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


pub fn prepare_energy_map<V>(index : V, energies : &Array1<f64>, unit : f64) -> FnvHashMap<i128, Vec<(V, usize)>>
    where V : EigenValue + Clone{
    // Prepare hashmap which will be used for degeneracy check

    let mut energy_map : FnvHashMap<i128, Vec<(V, usize)>> = FnvHashMap::default();
    for (idx, &e) in energies.iter().enumerate(){
        let k = (e / unit).floor() as i128;
        match energy_map.get_mut(&k){
            None => {
                match energy_map.get_mut(&(k -1)){
                    None => {
                        match energy_map.get_mut(&(k + 1)){
                            None => {
                                energy_map.insert(k, vec![(index, idx)]);
                            },
                            Some(v) => {
                                v.push((index, idx));
                            }
                        };
                    },
                    Some(v) => {
                        v.push((index, idx));
                    }
                };
            },
            Some(v) => {
                v.push((index, idx));
            },
        }
    }

    return energy_map;
}

pub fn count_degeneracy_from<V>(energy_map : &mut FnvHashMap<i128, Vec<(V, usize)>>, index : V, energies : &Array1<f64>, unit : f64) -> bool
    where V : EigenValue + Clone  + std::fmt::Debug{
    // Store degeneracy information only for value already in energy_map

    let mut t = false;
    for (idx, &e) in energies.iter().enumerate(){
        let k = (e / unit).floor() as i128;
        match energy_map.get_mut(&k){
            None => {
                match energy_map.get_mut(&(k - 1)){
                    None => {
                        match energy_map.get_mut(&(k + 1)){
                            None => {},
                            Some(v) => {
                                v.push((index, idx));
                                t = true;
                            }
                        };
                    },
                    Some(v) => {
                        v.push((index, idx));
                        t = true;
                    }
                };
            },
            Some(v) => {
                v.push((index, idx));
                t = true;
            },
        }
    }
    return t;
}

pub fn degeneracy_pair<V>(length : usize, energy_map : &FnvHashMap<i128, Vec<(V, usize)>>) -> Result<(Array1<f64>, FnvHashMap<V, Vec<(usize, usize)>>), Error>
    where V : EigenValue + Hash + Eq{
    let mut pair_map : FnvHashMap<V, Vec<(usize, usize)>> = FnvHashMap::default();
    let mut counts : Array1<f64> = Array1::ones(length);

    for (_i, vec) in energy_map.iter(){
        let (egn_v0, idx0) = vec[0];
        if length <= idx0{
            return Err(Error::make_error_syntax(ErrorCode::OverFlow));
        }

        counts[idx0] = 1.0 / (vec.len() as f64);
        for (egn_v, idx) in vec.iter().skip(1){
            if *egn_v == egn_v0{
                counts[*idx] = -1f64;
            }
            match pair_map.get_mut(egn_v){
                None => {
                    pair_map.insert(*egn_v, vec![(*idx, idx0)]);
                },
                Some(list) => {
                    list.push((*idx, idx0));
                }
            }
        }
    }

    return Ok((counts, pair_map));
}

pub fn degeneracy_triple<V>(length : usize, energy_map : &FnvHashMap<i128, Vec<(V, usize)>>) -> Result<(Array1<f64>, FnvHashMap<V, Vec<(usize, usize, usize)>>), Error>
    where V : EigenValue + Hash + Eq{
    let mut pair_map : FnvHashMap<V, Vec<(usize, usize, usize)>> = FnvHashMap::default();
    let mut counts : Array1<f64> = Array1::ones(length);

    for (_i, vec) in energy_map.iter(){
        let (egn_v0, idx0) = vec[0];
        if length <= idx0{
            return Err(Error::make_error_syntax(ErrorCode::OverFlow));
        }

        counts[idx0] = 1.0 / (vec.len() as f64);
        for (egn_v1, idx1) in vec.iter().skip(1){
            if *egn_v1 == egn_v0{
                counts[*idx1] = -1f64;
            }
            for (egn_v2, idx2) in vec.iter(){
                if *egn_v1 != *egn_v2{
                    continue;
                }
                match pair_map.get_mut(egn_v1){
                    None => {
                        pair_map.insert(*egn_v1, vec![(*idx1, *idx2, idx0)]);
                    },
                    Some(list) => {
                        list.push((*idx1, *idx2, idx0));
                    }
                }
            }
        }
    }

    return Ok((counts, pair_map));
}


#[cfg(test)]
mod test {
    use ndarray_linalg::close_l2;
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

    #[test]
    fn test_count_degeneracy_from(){
        let energies = arr1(&[1.11, 2.999, 3.0, 6.0, 5.99999999]);
        let mut energy_map = prepare_energy_map(0, &energies, 0.1);

        for (i, vector) in energy_map.iter(){
            match i{
                11 => {
                    assert_eq!(vector, &vec![(0, 0)]);
                },
                29 => {
                    assert_eq!(vector, &vec![(0, 1), (0, 2)]);
                },
                60 => {
                    assert_eq!(vector, &vec![(0, 3), (0, 4)]);
                },
                _ => panic!(),
            }
        }

        let others = arr1(&[1.0999, 2.0, 2.8, 2.91, 3.01, 3.11, 5.91, 6.01]);
        assert!(count_degeneracy_from(&mut energy_map, 1, &others, 0.1));

        for (i, vector) in energy_map.iter(){
            match i{
                11 => {
                    assert_eq!(vector, &vec![(0, 0), (1, 0)]);
                },
                29 => {
                    assert_eq!(vector, &vec![(0, 1), (0, 2), (1, 3), (1, 4)]);
                },
                60 => {
                    assert_eq!(vector, &vec![(0, 3), (0, 4), (1, 6), (1, 7)]);
                },
                _ => panic!(),
            }
        }
    }

    #[test]
    fn test_degeneracy_pair(){
        let length = 5;
        let mut energy_map : FnvHashMap<i128, Vec<(usize, usize)>> = FnvHashMap::default();
        energy_map.insert(0, vec![(1, 0), (1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 2), (3, 3)]);

        let (counts, pair_map) = degeneracy_pair(length, &energy_map).unwrap();
        close_l2(&counts, &arr1(&vec![0.111111111111111, -1.0, -1.0, -1.0, 1.0]), 1e-3);
        for (egn_v, pairs) in pair_map.iter(){
            match egn_v{
                1 => {
                    assert_eq!(pairs, &vec![(1, 0), (2, 0), (3, 0)]);
                },
                2 => {
                    assert_eq!(pairs, &vec![(1, 0), (2, 0), (3, 0)]);
                },
                3 => {
                    assert_eq!(pairs, &vec![(2, 0), (3, 0)]);
                },
                _ => panic!(),
            };
        }
    }

    #[test]
    fn test_degeneracy_triple(){
        let length = 5;
        let mut energy_map : FnvHashMap<i128, Vec<(usize, usize)>> = FnvHashMap::default();
        energy_map.insert(0, vec![(1, 0), (1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 2), (3, 3)]);

        let (counts, pair_map) = degeneracy_triple(length, &energy_map).unwrap();
        close_l2(&counts, &arr1(&vec![0.11111111111111, -1.0, -1.0, -1.0, 1.0]), 1e-3);
        for (egn_v, pairs) in pair_map.iter(){
            match egn_v{
                1 => {
                    assert_eq!(pairs, &vec![(1, 0, 0), (1, 1, 0), (1, 2, 0), (1, 3, 0),
                                           (2, 0, 0), (2, 1, 0), (2, 2, 0), (2, 3, 0),
                                           (3, 0, 0), (3, 1, 0), (3, 2, 0), (3, 3, 0)]);
                },
                2 => {
                    assert_eq!(pairs, &vec![(1, 1, 0), (1, 2, 0), (1, 3, 0),
                                           (2, 1, 0), (2, 2, 0), (2, 3, 0),
                                           (3, 1, 0), (3, 2, 0), (3, 3, 0)]);
                },
                3 => {
                    assert_eq!(pairs, &vec![(2, 2, 0), (2, 3, 0),
                                           (3, 2, 0), (3, 3, 0)]);
                },
                _ => panic!(),
            };
        }
    }
}



