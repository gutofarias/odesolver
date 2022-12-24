use odesolver::solver_vector as SV;
// use odesolver::solver_trait as ST;


// const ORDER:usize = 2;
// const ORDER_T:usize = ORDER + 1;

fn main() {
    let initial_state = vec!(3.0,2.0);
    // let initial_state = [3.0,2.0];
    let tstart = 0.0;
    let tend = 400.0;
    let step = 0.0001;
    let ratio_step_output = 100;
    let odeparam = SV::ODEParam {time : tstart, tend
                          ,step
                          ,ratio_step_output
                          };

    let file_vec = "./teste_vec.txt".to_string();
    let (data,_) = SV::solve_ode(system_function, odeparam, initial_state, SV::ODESolver::RK4); 
    
    SV::data_to_file(&data, file_vec, None).unwrap();

    
    // let sist = FSist1 { state : initial_state
    //                     // , outros_parametros : false
    // };

    // let file_trait = "./teste_trait.txt".to_string();
    // let(data,_,_) = ST::solve_ode::<ORDER,ORDER_T,_>(sist, odeparam, ST::ODESolver::RK4);
    // ST::data_to_file(&data, file_trait, None).unwrap();
    
}


// #[derive(Clone)]
// struct FSist1<const N: usize> {
//     state : ST::State<N>,
// }


fn system_function (_time: f64, state: &SV::State) -> SV::DState {
    let mut dstate =  vec!(0.0;2);
    dstate[0] = -0.5*state[0];
    dstate[1] = -0.00001*state[1];
    
    dstate
}

// impl ST::ODESystem<ORDER> for FSist1<ORDER> {
//     fn state (&self) -> &ST::State<ORDER>{
//         return &self.state;
//     }
    
//     fn dstate (&self, _time : f64) -> ST::DState<ORDER>{
//         let state = self.state;
//         let mut dstate  = [0.0; ORDER];
        
//         dstate[0] = -0.5*state[0];
//         dstate[1] =  -0.00001*state[1];

//         dstate
//     }

//     // fn new_from_state (&self, state : ST::State<ORDER>) -> Self {
//     //     FSist1::<ORDER> { state, .. *self}
//     // }

//     fn update_state(&mut self, state : ST::State<ORDER>) {
//         self.state = state;
//     }
// }
