pub use super::*;

/// Vector of the state of the system
pub type State = Vec<f64>;
/// Vector of the differential of the state of the system
pub type DState = Vec<f64>;
/// 2D vector of the data of the simulation
pub type Data = Vec<Vec<f64>>;

/// Function that receives as input the time and state of the system and returns the differential of the state.
pub type SystemFunction = fn (time: f64, state: &State) -> DState;
// pub type SystemFunction2<T: Fn(f64, &State) -> DState> = fn (func : T) -> DState;

// type Solver = fn (system_function : SystemFunction, step : f64, time : f64, state : &State) -> State;
type SolverClosure<T> = fn (system_function : &mut T, step : f64, time : f64, state : &State) -> State;

/// Solves one `step` of the ODE problem.
///
///# Inputs
///
///`system_function: Closure FnMut (f64,&State) -> DState` Also an closure Fn or an fn of the type `SystemFunction` are accepted. The input is a closure or fn which receives `time : f64` and `state : &State` returning a `DState`. In other words, function with the right parameters in order to obtain the differential of the state of the system.
///
///`odeparam: ODEParam`. A value of type ODEParam.
///
///`state: State`. A vector with an initial state.
///
///`odesolver: ODESolver`. A choice of an ODE solver.
pub fn edo_step<SysFunc : Fn (f64, &State) -> DState > (odeparam : ODEParam, system_function : &mut SysFunc, state : State, odesolver: ODESolver) -> State {
    
    let solver = match odesolver {
        ODESolver::RK4 => rk4_closure::<SysFunc>,
        ODESolver::Euler => euler_closure::<SysFunc>,
    };

    solver(system_function,odeparam.step
                    ,odeparam.time,&state)
}


// fn rk4 (system_function : SystemFunction, step : f64, time : f64, state : &State) -> State {
    
//     let k1 = system_function (time, state);
    
//     let time2 = time + 0.5*step;
//     let xs2 : Vec<f64> = state.iter().zip(&k1).map(|(x,k)|(*x) + 0.5 * k * step).collect();
//     let k2 = system_function(time2, &xs2);

//     let time3 = time2;
//     let xs3: Vec<f64> = state.iter().zip(&k2)
//         .map(|(x,kb)| (*x) + 0.5 * kb * step)
//         .collect();
//     let k3 = system_function(time3, &xs3);
    
//     let time4 = time + step;
//     let xs4:Vec<f64> = state.iter().zip(&k3)
//         .map(|(x,kc)| (*x) + kc * step).collect();
//     let k4:Vec<f64> = system_function(time4, &xs4); 
    
//     let mut output = Vec::<f64>::new(); 
//     for (i,_ka) in k1.iter().enumerate(){
//         output.push(state[i] + step*(1.0/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]));
//     }
    
//     output
// }

// fn euler (system_function : SystemFunction, step : f64, time : f64, state : &State) -> State {
//     let mut xps = system_function(time, &state);
//     for (x,xp) in state.iter().zip(xps.iter_mut()) {
//         *xp = (*x) + step*(*xp);
//         //tou aproveitando o proprio xp pra nao criar outra variavel pro valor final.
//     }

//     xps
// }


// pub type SystemFunction2<T: Fn(f64, &State) -> DState> = fn (func : T) -> DState;


fn rk4_closure<SysFunction : FnMut (f64,&State) -> DState> (system_function : &mut SysFunction, step : f64, time : f64, state : &State) -> State {
    
    let k1 = system_function (time, state);
    
    let time2 = time + 0.5*step;
    let xs2 : Vec<f64> = state.iter().zip(&k1).map(|(x,k)|(*x) + 0.5 * k * step).collect();
    let k2 = system_function(time2, &xs2);

    let time3 = time2;
    let xs3: Vec<f64> = state.iter().zip(&k2)
        .map(|(x,kb)| (*x) + 0.5 * kb * step)
        .collect();
    let k3 = system_function(time3, &xs3);
    
    let time4 = time + step;
    let xs4:Vec<f64> = state.iter().zip(&k3)
        .map(|(x,kc)| (*x) + kc * step).collect();
    let k4:Vec<f64> = system_function(time4, &xs4); 
    
    let mut output = Vec::<f64>::new(); 
    for (i,_ka) in k1.iter().enumerate(){
        output.push(state[i] + step*(1.0/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]));
    }
    
    output
}

fn euler_closure<SysFunction : FnMut (f64,&State) -> DState> (system_function : &mut SysFunction, step : f64, time : f64, state : &State) -> State {
    let mut xps = system_function(time, state);
    for (x,xp) in state.iter().zip(xps.iter_mut()) {
        *xp = (*x) + step*(*xp);
        //tou aproveitando o proprio xp pra nao criar outra variavel pro valor final.
    }

    xps
}



/// Main function for solving ODEs using vectors. Returns a tuple with (Data, ODEParam), where ODEParam have updated values.
///
///# Inputs
///
///`system_function: Closure FnMut (f64,&State) -> DState` Also an closure Fn or an fn of the type `SystemFunction` are accepted. The input is a closure or fn which receives `time : f64` and `state : &State` returning a `DState`. In other words, function with the right parameters in order to obtain the differential of the state of the system.
///
///`odeparam: ODEParam`. A value of type ODEParam.
///
///`odesolver: ODESolver`. A choice of an ODE solver.
pub fn solve_ode <SysFunc : FnMut (f64, &State) -> DState> (mut system_function: SysFunc, odeparam : ODEParam, state : State, odesolver: ODESolver) -> (Data, ODEParam) {

    let mut data : Vec<Vec<f64>> = Vec::<Vec<f64>>::new();
    
    let tini = odeparam.time;
    let tend = odeparam.tend;
    // let step = odeparam.step;
    // let ratio_step_output = odeparam.ratio_step_output;


    let mut datum: Vec<f64> = vec!(tini);
    datum.append(&mut state.clone());
    
    data.push(datum);

    let mut new_param : ODEParam = odeparam;
    let mut new_state : State = state;
    
    let solver = match odesolver {
        ODESolver::RK4 => rk4_closure,
        ODESolver::Euler => euler_closure,
    };
    
    loop {
        (new_state, new_param) =
            integrator_closure(new_param, &mut system_function, new_state, solver);

        let new_time = new_param.time;

        datum = vec!(new_time);
        datum.append(&mut new_state.clone());
        data.push(datum);

        if (new_time - tend).abs() < 1.0e-10 {
            break;
        }
    }

    (data, new_param)
}



fn integrator_closure<SysFunc : FnMut (f64, &State) -> DState> (odeparam : ODEParam, system_function : &mut SysFunc, state : State, solver : SolverClosure<SysFunc>) -> (Vec<f64>, ODEParam)
{
    let mut time = odeparam.time;
    let tend = odeparam.tend;
    let step = odeparam.step;
    let ratio_step_output = odeparam.ratio_step_output;
    // let solver = odeparam.solver;
    let mut new_state: Vec<f64> = state;

    for _i in  0 .. ratio_step_output {
        if (time + step) - tend > 0.0  {
            let new_step = tend - time;
            // println!("{}",new_step);
            new_state = (solver)(system_function, new_step, time, &new_state);
            time = tend;
            break;
        } else {
            new_state = (solver)(system_function, step, time, &new_state);
            time += step;
        }
    }

    let new_param = ODEParam {
        time, tend, step, ratio_step_output};
    
    (new_state, new_param)

}






/// Saves the data to a given filename/filepath
///
/// The first column of data is the times, the second onwards are the values of the state at that particular time. It has an `header: Option<String>` that when given a Some(String) will add the string as a header in the data file.
pub fn data_to_file(data:&Data,file_as_string:String, header : Option<String>) -> Result<(), Box<dyn Error>> {
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


/// Solves the ODE but at each iteration exports the values to a file. So it does not keep a data vector in memory.
///
///# Inputs
///
///`system_function: Closure FnMut (f64,&State) -> DState` Also an closure Fn or an fn of the type `SystemFunction` are accepted. The input is a closure or fn which receives `time : f64` and `state : &State` returning a `DState`. In other words, function with the right parameters in order to obtain the differential of the state of the system.
///
///`odeparam: ODEParam`. A value of type ODEParam.
///
///`odesolver: ODESolver`. A choice of an ODE solver.
///
///`filestr: String`. String with a given filename/filepath to save the data
pub fn solve_ode_to_file <SysFunc : FnMut (f64,&State) -> DState> (mut system_function: SysFunc, odeparam : ODEParam, state : State, odesolver : ODESolver, filestr : String ) -> () {

    
    let tini = odeparam.time;
    let tend = odeparam.tend;
    // let step = odeparam.step;
    // let ratio_step_output = odeparam.ratio_step_output;


    let mut file = File::create(filestr).unwrap();
    
    let mut datum: Vec<f64> = vec!(tini);
    datum.append(&mut state.clone());
    
    
    let mut string = String::new();
    for value in datum {
            string += &format!("{:10.8e}\t",value).to_string();
    }
    writeln!(file,"{}", string).unwrap();

    let mut new_param : ODEParam = odeparam;
    let mut new_state : State = state;
    
    let solver = match odesolver {
        ODESolver::RK4 => rk4_closure,
        ODESolver::Euler => euler_closure,
    };
    
    loop {
        (new_state, new_param) =
            integrator_closure(new_param, &mut system_function, new_state, solver);

        let new_time = new_param.time;

        datum = vec!(new_time);
        datum.append(&mut new_state.clone());
        
        string.clear();
        for value in datum {
                string += &format!("{:10.8e}\t",value).to_string();
        }
        writeln!(file,"{}", string).unwrap();

        if (new_time - tend).abs() < 1.0e-10 {
            break;
        }
    }

}

