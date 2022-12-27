//! #ODESolver
//!
//! `odesolver` is a library for solving ODEs having two main ways of using it. One is by the module `solver_trait` which is faster but the functions need to know the number of variables at compile time. The other module is `solver_vector` it is not as fast as the other but the number of variables can be specified at run time.


use std::error::Error;
use std::fs::File;
use std::io::Write;



/// Struct ODEParam
///
/// It has all the parameters needed for the ODE solver: the actual `time`, the time for ending the simulation `tend`, the integration `step` and `relStepOut` which gives a ratio between the step time and the output time. If `ratio_step_output` is 1, at every step, output data will be added to the container, if it is 3, at every 3 steps, output data will be added to the container.
#[derive(Clone)]
pub struct ODEParam {
    pub time : f64,
    pub tend : f64,
    pub step : f64,
    pub ratio_step_output: usize,
}

/// Enum ODESolver
///
/// Contains the available solvers in the library
pub enum ODESolver {
    RK4,
    Euler,
}

/// odesolver using traits and known size at compile time
///
///# Example:
///
///```
///use odesolver::solver_trait as ST;
///
///const ORDER:usize = 2;
///const ORDER_T:usize = ORDER + 1;
/// // it is always ORDER + 1, so that the exporting Data has an extra field for the time.
///
///fn main() {
///    let initial_state = [3.0,2.0];
///    let tstart = 0.0;
///    let tend = 400.0;
///    let step = 0.0001;
///    let ratio_step_output = 100;
///    let odeparam = ST::ODEParam {time : tstart, tend
///                          ,step
///                          ,ratio_step_output
///                          };
///    
///    let sist = FSist1 { state : initial_state };
///
///    let file_trait = "./test_trait.txt".to_string();
///    let(data,_,_) = ST::solve_ode::<ORDER,ORDER_T,_>(sist, odeparam, ST::ODESolver::RK4);
///    
///    ST::data_to_file(&data, file_trait, None).unwrap();
///}
///
///
///#[derive(Clone)]
///struct FSist1<const N: usize> {
///    state : ST::State<N>,
///}
///
///impl ST::ODESystem<ORDER> for FSist1<ORDER> {
///    fn state (&self) -> &ST::State<ORDER>{
///        return &self.state;
///    }
///    
///    fn dstate (&self, _time : f64) -> ST::DState<ORDER>{
///        let state = self.state;
///        let mut dstate  = [0.0; ORDER];
///        
///        dstate[0] = -0.5*state[0];
///        dstate[1] =  -0.00001*state[1];
///
///        dstate
///    }
///
///    fn update_state(&mut self, state : ST::State<ORDER>) {
///        self.state = state;
///    }
///}
///```
pub mod solver_trait;






///odesolver using vectors, no need to know the size at compile time, but not as fast due to this.
///
///# Example
///
///```
///use odesolver::solver_vector as SV;
///
///fn main() {
///    let initial_state = vec!(3.0,2.0);
///    let tstart = 0.0;
///    let tend = 400.0;
///    let step = 0.0001;
///    let ratio_step_output = 100;
///    let odeparam = SV::ODEParam {time : tstart, tend
///                          ,step
///                          ,ratio_step_output
///                          };
///
///    let file_vec = "./test_vec.txt".to_string();
///    let (data,_) = SV::solve_ode(system_function, odeparam, initial_state, SV::ODESolver::RK4); 
///    
///    SV::data_to_file(&data, file_vec, None).unwrap();
///}
///
///fn system_function (_time: f64, state: &SV::State) -> SV::DState {
///    let mut dstate =  vec!(0.0;2);
///    dstate[0] = -0.5*state[0];
///    dstate[1] = -0.00001*state[1];
///    
///    dstate
///}
///```
pub mod solver_vector;






///odesolver using vectors but with trait bounds.
///
///It can be usefull if you are working with your own data types that happens to be defined by ODEs. In that case you just need to make your data type adhere to the trait in order to solve the ODEs, no matter how your data is organized. Also in this case there is no need to know the size at compile time, but not as fast due to this.
///
///# Example
///
///```
///use odesolver::solver_vector_trait as SVT;
///
///fn main() {
///    let initial_state = vec!(3.0,2.0);
///    // let initial_state = [3.0,2.0];
///    let tstart = 0.0;
///    let tend = 400.0;
///    let step = 0.0001;
///    let ratio_step_output = 100;
///    let odeparam = SVT::ODEParam {time : tstart, tend
///                          ,step
///                          ,ratio_step_output
///                          };
///
///    let sist = FSist2 { state : initial_state };
///    let file_vec_trait = "./teste_vec_trait.txt".to_string();
///    let(data,_,_) = SVT::solve_ode::<_>(sist, odeparam, SVT::ODESolver::RK4);
///    
///    SVT::data_to_file(&data, file_vec_trait, None).unwrap();
///}
///
///#[derive(Clone)]
///struct FSist2 {
///    state : SVT::State,
///}
///
///impl SVT::ODESystem for FSist2 {
///    fn state (&self) -> &SVT::State{
///        return &self.state;
///    }
///    
///    fn dstate (&self, _time : f64) -> SVT::DState{
///        let state = self.state();
///        let mut dstate  = Vec::<f64>::new();
///        
///        dstate.push(-0.5*state[0]); //[0]
///        dstate.push(-0.00001*state[1]); //[1]
///
///        dstate
///    }
///
///    fn update_state(&mut self, state : SVT::State) {
///        self.state = state;
///    }
///}
///```
pub mod solver_vector_trait;
