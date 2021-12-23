use crate::bases::BasisGenerator;
use crate::states::{SimpleState, momentum::{EigenNumMomentum, NumMomentumState}};
use fnv::FnvHashMap;

pub type BasisNK = BasisGenerator<EigenNumMomentum>;

impl BasisNK{
    pub fn new(v : EigenNumMomentum, length : usize) -> Self{
        Self{
            length,
            value : Box::new(v),
        }
    }

    pub fn check_commensurability(&self, period : usize) -> bool{
        self.value.check_commensurability(period, self.length)
    }

    pub fn build(&self) -> (Vec<NumMomentumState>, FnvHashMap<usize, (usize, usize)>){
        let max_state = 1 << self.length;
        let mut basis : Vec<NumMomentumState> = Vec::new();
        let mut indices : FnvHashMap<usize, (usize, usize)> = FnvHashMap::default();
        let mut idx = 0;

        for n in 0..max_state{
            let state = SimpleState::new(n, self.length);
            let m = state.total_number();
            let p = state.period();

            if m != self.value.total_number() {
                continue;
            }

            if let Some(nkstate) = NumMomentumState::new(state, self.value.wave_number()){
                for (i, state) in nkstate.state.states.iter().enumerate(){
                    let num = state.rep;
                    indices.insert(num, (idx, (p - i) % p));
                }

                basis.push(nkstate);
                idx += 1;
            }
        }

        return (basis, indices);
    }
}

