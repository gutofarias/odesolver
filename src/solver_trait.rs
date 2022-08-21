pub use super::*;

pub type State<const N: usize> = [f64;N];
pub type DState<const N: usize> = [f64;N];
pub type Data<const M: usize> = Vec<[f64;M]>;

pub trait EdoSystem<const N: usize> {
    
    fn state (&self) -> &State<N>; 

    fn dstate (&self, time : f64) -> DState<N>;

    fn new_from_state (&self, state : State<N>) -> Self; 

    // fn update_state (&mut self, state : State<N>);
}


pub type Solver<const N:usize, Sist: EdoSystem<N>> = fn (sist: &Sist, step : f64, time : f64) -> State<N>;


fn integrator<const N:usize, Sist: EdoSystem<N>> (sist : Sist, edoParam : EdoParam, solver : Solver<N,Sist>) -> (Sist, EdoParam)
{
    let mut time = edoParam.time;
    let tend = edoParam.tend;
    let step = edoPakam.step;
    let relStepOut = edoParam.relStepOut;
    // let solver = edoParam.solver;
    let mut newState = sist.state().clone();
    let mut newSist = sist; 

    for i in  0 .. relStepOut {
        if (time + step) - tend > 0.0  {
            let new_step = tend - time;
            // println!("{}",new_step);
            newState = (solver)(&newSist, step, time);
            newSist = newSist.new_from_state(newState);
            time = tend;
            break;
        } else {
            newState = (solver)(&newSist, step, time);
            newSist = newSist.new_from_state(newState);
            time += step;
        }
    }

    let newParam = EdoParam {
        time, .. edoParam};
    
    (newSist, newParam)
}



pub fn edoSolver<const N:usize, const M:usize, Sist: EdoSystem<N>> (sist : Sist, edoparam : EdoParam, solver : Solver<N,Sist>) -> (Data<M>, Sist, EdoParam) {
    
    let tini = edoparam.time;
    let tend = edoparam.tend;
    let step = edoparam.step;
    let relStepOut = edoparam.relStepOut;

    let capacity = ((tend-tini)/(step*(relStepOut as f64)) + 10.0) as usize;
    
    let mut data : Vec<[f64;M]> = Vec::<[f64;M]>::with_capacity(capacity);

    let mut datum: [f64;M] = [0.0;M];

    datum[0] = tini;
    for (i,x) in sist.state().clone().iter().enumerate() {
        datum[i+1] = *x;
    }
                                                  
        // = vec!(tini);
    // datum.append(&mut state.clone());
    
    data.push(datum);

    let mut newParam  = edoparam;
    let mut newSist  = sist;
    loop {
        (newSist, newParam) =
            integrator::<N,Sist>(newSist,newParam,solver);

        let newTime = newParam.time;

        // datum = vec!(newTime);
        // datum.append(&mut newState.clone());
        
        datum[0] = newTime;
        for (i,x) in newSist.state().clone().iter().enumerate() {
            datum[i+1] = *x;
        }
        
        data.push(datum);

        if (newTime - tend).abs() < 1.0e-10 {
            break;
        }
    }

    (data, newSist, newParam)
}




pub fn euler<const N:usize, Sist: EdoSystem<N>> (sist: &Sist, step : f64, time : f64) -> State<N> {
    let state = sist.state();
    let mut xps = sist.dstate(time);// funcSist(time, &state);
    for (x,xp) in state.iter().zip(xps.iter_mut()) {
        *xp = (*x) + step*(*xp);
        //tou aproveitando o proprio xp pra nao criar outra variavel pro valor final.
    }

    xps
}



pub fn rk4<const N:usize, Sist: EdoSystem<N>> (sist: &Sist, step : f64, time : f64) -> State<N> {    
    
    let state = sist.state();
    
    let k1 = sist.dstate(time);

    let time2 = time + 0.5*step;
    let mut xs2  = [0.0;N];
    for (i,x) in state.iter().enumerate(){
        xs2[i] = *x + 0.5*k1[i]*step;
    }
    let sist2 = sist.new_from_state(xs2);
    let k2 = sist2.dstate(time2);

    let time3 = time2;
    let mut xs3 = [0.0;N];
    for (i,x) in state.iter().enumerate(){
        xs3[i] = *x + 0.5*k2[i]*step;
    }
    let sist3 = sist.new_from_state(xs3);
    let k3 = sist3.dstate(time3);
    
    let time4 = time + step;
    let mut xs4 = [0.0;N];
    for (i,x) in state.iter().enumerate(){
        xs4[i] = *x + k3[i]*step;
    }
    let sist4 = sist.new_from_state(xs4);
    let k4 = sist4.dstate(time4); 
    
    let state = sist.state();
    let mut output = [-15.0; N]; 
    for (i,_ka) in k1.iter().enumerate(){
        output[i] = state[i] + step*(1.0/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
    
    output
}


// use super::solver_const::{data_to_file_const, DataConst};

// pub fn data_to_file_trait(data:&DataConst,file_as_string:String) -> Result<(), Box<dyn Error>> {
//     data_to_file_const(data,file_as_string)
// }

pub fn data_to_file<const M : usize>(data:&Data<M>,file_as_string:String, header : Option<String>) -> Result<(), Box<dyn Error>> {
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

