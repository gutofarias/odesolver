pub use super::*;

pub type State = Vec<f64>;
pub type DState = Vec<f64>;
pub type Data = Vec<Vec<f64>>;

pub type FuncSist = fn (time: f64, state: &State) -> DState;

pub type Solver = fn (funcSist : FuncSist, step : f64, time : f64, state : &State) -> State;


pub fn edoStep (edoparam : EdoParam, funcSist : FuncSist, state : State, solver: Solver) -> State {

    solver(funcSist,edoparam.step
                    ,edoparam.time,&state)
}


pub fn rk4 (funcSist : FuncSist, step : f64, time : f64, state : &State) -> State {
    
    let k1 = funcSist (time, state);
    
    let time2 = time + 0.5*step;
    let xs2 : Vec<f64> = state.iter().zip(&k1).map(|(x,k)|(*x) + 0.5 * k * step).collect();
    let k2 = funcSist(time2, &xs2);

    let time3 = time2;
    let xs3: Vec<f64> = state.iter().zip(&k2)
        .map(|(x,kb)| (*x) + 0.5 * kb * step)
        .collect();
    let k3 = funcSist(time3, &xs3);
    
    let time4 = time + step;
    let xs4:Vec<f64> = state.iter().zip(&k3)
        .map(|(x,kc)| (*x) + kc * step).collect();
    let k4:Vec<f64> = funcSist(time4, &xs4); 
    
    let mut output = Vec::<f64>::new(); 
    for (i,_ka) in k1.iter().enumerate(){
        output.push(state[i] + step*(1.0/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]));
    }
    
    output
}

pub fn euler (funcSist : FuncSist, step : f64, time : f64, state : &State) -> State {
    let mut xps = funcSist(time, &state);
    for (x,xp) in state.iter().zip(xps.iter_mut()) {
        *xp = (*x) + step*(*xp);
        //tou aproveitando o proprio xp pra nao criar outra variavel pro valor final.
    }

    xps
}




pub fn edoSolver (funcSist: FuncSist, edoparam : EdoParam, state : State, solver: Solver) -> (Data, EdoParam) {

    let mut data : Vec<Vec<f64>> = Vec::<Vec<f64>>::new();
    
    let tini = edoparam.time;
    let tend = edoparam.tend;
    let step = edoparam.step;
    let relStepOut = edoparam.relStepOut;


    let mut datum: Vec<f64> = vec!(tini);
    datum.append(&mut state.clone());
    
    data.push(datum);

    let mut newParam : EdoParam = edoparam;
    let mut newState : State = state;
    loop {
        (newState, newParam) =
            integrator(newParam, funcSist, newState, solver);

        let newTime = newParam.time;

        datum = vec!(newTime);
        datum.append(&mut newState.clone());
        data.push(datum);

        if (newTime - tend).abs() < 1.0e-10 {
            break;
        }
    }

    (data, newParam)
}




fn integrator ( edoParam : EdoParam, funcSist : FuncSist, state : State, solver : Solver) -> (Vec<f64>, EdoParam)
{
    let mut time = edoParam.time;
    let tend = edoParam.tend;
    let step = edoParam.step;
    let relStepOut = edoParam.relStepOut;
    // let solver = edoParam.solver;
    let mut newState: Vec<f64> = state;

    for i in  0 .. relStepOut {
        if (time + step) - tend > 0.0  {
            let new_step = tend - time;
            // println!("{}",new_step);
            newState = (solver)(funcSist, new_step, time, &newState);
            time = tend;
            break;
        } else {
            newState = (solver)(funcSist, step, time, &newState);
            time += step;
        }
    }

    let newParam = EdoParam {
        time, tend, step, relStepOut};
    
    (newState, newParam)

}



pub fn data_to_file(data:&Data,file_as_string:String, header : Option<String>) -> Result<(), Box<dyn Error>> {
    let file = File::create(file_as_string);

    let mut file = match file {
        Ok(f) => f,
        Err(erro) => {
            eprintln!("Não foi possível criar o arquivo para exportar os dados. Erro: {}", erro);
            return Err(Box::new(erro) as Box<dyn Error>);},
    };
    
    if let Some(headerstring) = header {
        match writeln!(file,"{}", headerstring){
            Ok(_) => (),
            Err(erro) => {
            eprintln!("Erro ao escrever no arquivo. Erro: {}", erro);
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
            eprintln!("Erro ao escrever no arquivo. Erro: {}", erro);
            return Err(Box::new(erro) as Box<dyn Error>);},
        }
    }

    Ok(())
}














pub fn edoSolverToFile (funcSist: FuncSist, edoparam : EdoParam, state : State, filestr : String, solver: Solver) -> () {

    
    let tini = edoparam.time;
    let tend = edoparam.tend;
    let step = edoparam.step;
    let relStepOut = edoparam.relStepOut;


    let mut file = File::create(filestr).unwrap();
    
    let mut datum: Vec<f64> = vec!(tini);
    datum.append(&mut state.clone());
    
    
    let mut string = String::new();
    for value in datum {
            string += &format!("{:10.8e}\t",value).to_string();
    }
    writeln!(file,"{}", string).unwrap();

    let mut newParam : EdoParam = edoparam;
    let mut newState : State = state;
    loop {
        (newState, newParam) =
            integrator(newParam, funcSist, newState, solver);

        let newTime = newParam.time;

        datum = vec!(newTime);
        datum.append(&mut newState.clone());
        
        string.clear();
        for value in datum {
                string += &format!("{:10.8e}\t",value).to_string();
        }
        writeln!(file,"{}", string).unwrap();

        if (newTime - tend).abs() < 1.0e-10 {
            break;
        }
    }

}

