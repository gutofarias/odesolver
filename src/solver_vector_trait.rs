pub use super::*;

/// Vector of the state of the system
pub type State = Vec<f64>;
/// Vector of the differential of the state of the system
pub type DState = Vec<f64>;
/// 2D vector of the data of the simulation
pub type Data = Vec<Vec<f64>>;


///Trait used to caracterize a data type as being a system defined by ODEs. 
///
///If this trait is implemented you can use it to solve the ODEs on your own data type.
pub trait ODESystem {
    
    ///Return the actual state of the system.
    fn state (&self) -> &State; 

    ///Returns the differential of the state of the system.
    fn dstate (&self, time : f64) -> DState;

    ///Updates the state of the system.
    fn update_state (&mut self, state : State);
}


type Solver<Sist> = fn (sist: &Sist, step : f64, time : f64) -> State;


fn integrator<Sist: ODESystem> (sist : Sist, odeparam : ODEParam, solver : Solver<Sist>) -> (Sist, ODEParam)
{
    let mut time = odeparam.time;
    let tend = odeparam.tend;
    let step = odeparam.step;
    let ratio_step_output = odeparam.ratio_step_output;
    // let solver = odeparam.solver;
    let mut new_state;  // = sist.state().clone();
    let mut new_sist = sist; 

    for _i in  0 .. ratio_step_output {
        if (time + step) - tend > 0.0  {
            let new_step = tend - time;
            new_state = (solver)(&new_sist, new_step, time);
            new_sist.update_state(new_state);
            time = tend;
            break;
        } else {
            new_state = (solver)(&new_sist, step, time);
            new_sist.update_state(new_state);
            time += step;
        }
    }

    let new_param = ODEParam {
        time, .. odeparam};
    
    (new_sist, new_param)
}

              
/// Main function for solving ODEs, returns an tuple with (Data, Sist, ODEParam) where Sist and ODEParam have updated values.
///
///# inputs
///
///`sist : Sist`. Any type which implements ODESystem and Clone
///
///`odeparam : ODEParam`. An value of ODEParam.
///
///`odesolver : ODESolver`. A choice of an ODE solver.
pub fn solve_ode<Sist: ODESystem + Clone> (sist : Sist, odeparam : ODEParam, odesolver: ODESolver) -> (Data, Sist, ODEParam) {
    
    let tini = odeparam.time;
    let tend = odeparam.tend;
    let step = odeparam.step;
    let ratio_step_output = odeparam.ratio_step_output;

    let capacity = ((tend-tini)/(step*(ratio_step_output as f64)) + 10.0) as usize;
    
    let mut data : Vec<Vec<f64>> = Vec::<Vec<f64>>::with_capacity(capacity);

    let mut datum: Vec<f64> = vec!(tini);
    datum.append(&mut sist.state().clone());
    
    data.push(datum);
    
    let mut new_param  = odeparam;
    let mut new_sist  = sist;
    
    let solver = match odesolver {
        ODESolver::RK4 => rk4,
        ODESolver::Euler => euler,
    };
    
    loop {
        (new_sist, new_param) =
            integrator::<Sist>(new_sist,new_param,solver);

        let new_time = new_param.time;

        datum = vec!(new_time);
        datum.append(&mut new_sist.state().clone());
        
        data.push(datum);

        if (new_time - tend).abs() < 1.0e-10 {
            break;
        }
    }

    (data, new_sist, new_param)
}



fn euler<Sist: ODESystem> (sist: &Sist, step : f64, time : f64) -> State {
    let state = sist.state();
    let mut xps = sist.dstate(time);// funcSist(time, &state);
    for (x,xp) in state.iter().zip(xps.iter_mut()) {
        *xp = (*x) + step*(*xp);
    }

    xps
}


fn rk4<Sist: ODESystem + Clone> (sist: &Sist, step : f64, time : f64) -> State {    
    
    let state = sist.state();
    let k1 = sist.dstate(time);
    
    
    let time2 = time + 0.5*step;
    let xs2 : Vec<f64> = state.iter().zip(&k1).map(|(x,k)|(*x) + 0.5 * k * step).collect();

    let mut sist2 = sist.clone();
    sist2.update_state(xs2);
    let k2 = sist2.dstate(time2);

    let time3 = time2;
    let xs3: Vec<f64> = state.iter().zip(&k2)
        .map(|(x,kb)| (*x) + 0.5 * kb * step)
        .collect();
    
    let mut sist3 = sist.clone();
    sist3.update_state(xs3);
    let k3 = sist3.dstate(time3);
    
    let time4 = time + step;
    let xs4:Vec<f64> = state.iter().zip(&k3)
        .map(|(x,kc)| (*x) + kc * step).collect();
    
    let mut sist4 = sist.clone();
    sist4.update_state(xs4);
    let k4 = sist4.dstate(time4); 
    
    let mut output = Vec::<f64>::new(); 
    for (i,_ka) in k1.iter().enumerate(){
        output.push(state[i] + step*(1.0/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]));
    }

    output
}


pub use crate::solver_vector::data_to_file;
