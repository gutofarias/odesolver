pub use super::*;

///Array of the state of a system
pub type State<const N: usize> = [f64;N];
///Array of the differential of the state of a system
pub type DState<const N: usize> = [f64;N];
///Array of data of the simulation
///
/// It should have size M = N+1. Rust does not allow for arithmetics on const generics, so I had to use an other generic.
pub type Data<const M: usize> = Vec<[f64;M]>;

///Trait used to caracterize a data type as being a system defined by ODEs. An order N (number of ODEs) has to be specified.
///
///If this trait is implemented you can use it to solve the ODEs on your own data type.
pub trait ODESystem<const N: usize> {
    
    ///Return the actual state of the system.
    fn state (&self) -> &State<N>; 

    ///Returns the differential of the state of the system.
    fn dstate (&self, time : f64) -> DState<N>;

    // fn new_from_state (&self, state : State<N>) -> Self; 

    ///Updates the state of the system.
    fn update_state (&mut self, state : State<N>);
}


// pub type Solver<const N:usize, Sist: EdoSystem<N>> = fn (sist: &Sist, step : f64, time : f64) -> State<N>;
// ODE Solver function
type Solver<const N:usize, Sist> = fn (sist: &Sist, step : f64, time : f64) -> State<N>;


fn integrator<const N:usize, Sist: ODESystem<N>> (sist : Sist, odeparam : ODEParam, solver : Solver<N,Sist>) -> (Sist, ODEParam)
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
            // new_state = (solver)(&new_sist, step, time);
            new_state = (solver)(&new_sist, new_step, time);
            // new_sist = new_sist.new_from_state(new_state);
            new_sist.update_state(new_state);
            time = tend;
            break;
        } else {
            new_state = (solver)(&new_sist, step, time);
            // new_sist = new_sist.new_from_state(new_state);
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
///`sist : Sist`. Any type which implements `ODESystem<N>` and `Clone`
///
///`odeparam : ODEParam`. An value of ODEParam.
///
///`odesolver : ODESolver`. A choice of an ODE solver.
pub fn solve_ode<const N:usize, const M:usize, Sist: ODESystem<N> + Clone> (sist : Sist, odeparam : ODEParam, odesolver: ODESolver) -> (Data<M>, Sist, ODEParam) {
    
    let tini = odeparam.time;
    let tend = odeparam.tend;
    let step = odeparam.step;
    let ratio_step_output = odeparam.ratio_step_output;

    let capacity = ((tend-tini)/(step*(ratio_step_output as f64)) + 10.0) as usize;
    
    let mut data : Vec<[f64;M]> = Vec::<[f64;M]>::with_capacity(capacity);

    let mut datum: [f64;M] = [0.0;M];

    datum[0] = tini;
    for (i,x) in sist.state().clone().iter().enumerate() {
        datum[i+1] = *x;
    }
                                                  
        // = vec!(tini);
    // datum.append(&mut state.clone());
    
    data.push(datum);

    let mut new_param  = odeparam;
    let mut new_sist  = sist;
    
    let solver = match odesolver {
        ODESolver::RK4 => rk4,
        ODESolver::Euler => euler,
    };
    
    loop {
        (new_sist, new_param) =
            integrator::<N,Sist>(new_sist,new_param,solver);

        let new_time = new_param.time;

        // datum = vec!(newTime);
        // datum.append(&mut newState.clone());
        
        datum[0] = new_time;
        for (i,x) in new_sist.state().clone().iter().enumerate() {
            datum[i+1] = *x;
        }
        
        data.push(datum);

        if (new_time - tend).abs() < 1.0e-10 {
            break;
        }
    }

    (data, new_sist, new_param)
}




fn euler<const N:usize, Sist: ODESystem<N>> (sist: &Sist, step : f64, time : f64) -> State<N> {
    let state = sist.state();
    let mut xps = sist.dstate(time);// funcSist(time, &state);
    for (x,xp) in state.iter().zip(xps.iter_mut()) {
        *xp = (*x) + step*(*xp);
        //tou aproveitando o proprio xp pra nao criar outra variavel pro valor final.
    }
    xps
}



fn rk4<const N:usize, Sist: ODESystem<N> + Clone> (sist: &Sist, step : f64, time : f64) -> State<N> {    
    
    let state = sist.state();
    
    let k1 = sist.dstate(time);

    let time2 = time + 0.5*step;
    let mut xs2  = [0.0;N];
    for (i,x) in state.iter().enumerate(){
        xs2[i] = *x + 0.5*k1[i]*step;
    }
    // let sist2 = sist.new_from_state(xs2);
    let mut sist2 = sist.clone();
    sist2.update_state(xs2);
    let k2 = sist2.dstate(time2);

    let time3 = time2;
    let mut xs3 = [0.0;N];
    for (i,x) in state.iter().enumerate(){
        xs3[i] = *x + 0.5*k2[i]*step;
    }
    // let sist3 = sist.new_from_state(xs3);
    let mut sist3 = sist.clone();
    sist3.update_state(xs3);
    let k3 = sist3.dstate(time3);
    
    let time4 = time + step;
    let mut xs4 = [0.0;N];
    for (i,x) in state.iter().enumerate(){
        xs4[i] = *x + k3[i]*step;
    }
    // let sist4 = sist.new_from_state(xs4);
    let mut sist4 = sist.clone();
    sist4.update_state(xs4);
    let k4 = sist4.dstate(time4); 
    
    let state = sist.state();
    let mut output = [0.0; N]; 
    for (i,_ka) in k1.iter().enumerate(){
        output[i] = state[i] + step*(1.0/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
    
    output
}


/// Saves the data to a given filename/filepath
///
/// The first column of data is the times, the second onwards are the values of the state at that particular time. It has an `header: Option<String>` that when given a Some(String) will add the string as a header in the data file.
pub fn data_to_file<const M : usize>(data:&Data<M>,file_as_string:String, header : Option<String>) -> Result<(), Box<dyn Error>> {
    let file = File::create(file_as_string);

    let mut file = match file {
        Ok(f) => f,
        Err(erro) => {
            eprintln!("Not possible to create/find the file for exporting the data. Error: {}", erro);
            return Err(Box::new(erro) as Box<dyn Error>);},
    };

    if let Some(headerstring) = header {
        match writeln!(file,"{}", headerstring){
            Ok(_) => (),
            Err(erro) => {
            eprintln!("Error in writing to file. Error: {}", erro);
            return Err(Box::new(erro) as Box<dyn Error>);},
        }
    }

    for datum in data {
        let mut string = String::new();
        for value in datum {
             string += &format!("{:10.8e}\t",value).to_string();
        }
        match writeln!(file,"{}", string) {
            Ok(_) => (),
            Err(erro) => {
            eprintln!("Error in writing to file. Error: {}", erro);
            return Err(Box::new(erro) as Box<dyn Error>);},
        }
    }

    Ok(())
}

