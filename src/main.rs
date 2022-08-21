// use rk4::*;
use ODESolver::solver_vector as SV;
use ODESolver::solver_trait as ST;


const order:usize = 2;
const order_t:usize = order + 1;

fn main() {
    let xs = vec!(3.0,2.0);
    let xsconst = [3.0,2.0];
    let tini = 0.0;
    let tfim = 400.0;
    let passo = 0.00001;
    let relSaida = 100;
    let param = ST::EdoParam {time : tini, tend : tfim
                          ,step : passo
                          ,relStepOut : relSaida
                          };
    let param2 = param.clone();

    let file_vec = "./teste_vec.txt".to_string();
    let (data,_) = SV::edoSolver(funcsist1, param, xs, SV::rk4); 
    SV::data_to_file(&data, file_vec, None).unwrap();

    
    let sist = FSist1 { state : [3.0,2.0]
                        , outros_parametros : false
    };

    let file_trait = "./teste_trait.txt".to_string();
    let(data,_,_) = ST::edoSolver::<order,order_t,_>(sist, param2, ST::rk4);
    ST::data_to_file(&data, file_trait, None).unwrap();

    
    
}

struct FSist1<const N: usize> {
    state : ST::State<N>,
    outros_parametros : bool,
}


fn funcsist1 (time: f64, state: &SV::State) -> SV::DState {
    let mut dstate =  vec!(0.0;2);
    dstate[0] = -0.5*state[0];
    dstate[1] = -0.00001*state[1];
    
    dstate
}

impl ST::EdoSystem<order> for FSist1<order> {
    fn state (&self) -> &ST::State<order>{
        return &self.state;
    }
    
    fn dstate (&self, time : f64) -> ST::DState<order>{
        let state = self.state;
        let mut dstate  = [0.0; order];
        
        dstate[0] = -0.5*state[0];
        dstate[1] =  -0.00001*state[1];

        dstate
    }

    fn new_from_state (&self, state : ST::State<order>) -> Self {
        FSist1::<order> { state, .. *self}
    }
}
