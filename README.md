**odesolver** is a **rust** library for solving ODE's.

The library has two possible interfaces, the first is module `solver_trait` in which the user knows the order of the problem (number of ODE's in a system of ODE's) at compile time. This interface uses rust arrays as data collection and traits to define the behavior of the system. As the compiler knows the sizes, it can allocate the needed memory in the stack, running the program faster. This module solves the problem roughly one order of magnitude faster (around 10 times) than the other module. 

The second is module `solver_vector` in which the user uses vectors as data collection and the memory is allocated at the heap. This allows greater flexibility for solving the problem as the order of the problem does not need to be known at compile time. 

# Examples

## Example `solver_trait`

```rust
use odesolver::solver_trait as ST;

const ORDER:usize = 2;
const ORDER_T:usize = ORDER + 1;
 //always ORDER + 1, so that the exporting Data has an extra field for the time.

fn main() {
    let initial_state = [3.0,2.0];
    let tstart = 0.0;
    let tend = 400.0;
    let step = 0.0001;
    let ratio_step_output = 100;
    let odeparam = ST::ODEParam {time : tstart, tend
                          ,step
                          ,ratio_step_output
                          };
    
    let sist = FSist1 { state : initial_state };

    let file_trait = "./test_trait.txt".to_string();
    let(data,_,_) = ST::solve_ode::<ORDER,ORDER_T,_>(sist, odeparam, ST::ODESolver::RK4);
    
    ST::data_to_file(&data, file_trait, None).unwrap();
}


#[derive(Clone)]
struct FSist1<const N: usize> {
    state : ST::State<N>,
}

impl ST::ODESystem<ORDER> for FSist1<ORDER> {
    fn state (&self) -> &ST::State<ORDER>{
        return &self.state;
    }
    
    fn dstate (&self, _time : f64) -> ST::DState<ORDER>{
        let state = self.state;
        let mut dstate  = [0.0; ORDER];
        
        dstate[0] = -0.5*state[0];
        dstate[1] =  -0.00001*state[1];

        dstate
    }

    fn update_state(&mut self, state : ST::State<ORDER>) {
        self.state = state;
    }
}
```

## Example `solver_vector`

```rust
use odesolver::solver_vector as SV;

fn main() {
    let initial_state = vec!(3.0,2.0);
    let tstart = 0.0;
    let tend = 400.0;
    let step = 0.0001;
    let ratio_step_output = 100;
    let odeparam = SV::ODEParam {time : tstart, tend
                          ,step
                          ,ratio_step_output
                          };

    let file_vec = "./test_vec.txt".to_string();
     (data,_) = SV::solve_ode(system_function, odeparam, initial_state, SV::ODESolver::RK4); 
    
    SV::data_to_file(&data, file_vec, None).unwrap();
}

fn system_function (_time: f64, state: &SV::State) -> SV::DState {
    let mut dstate =  vec!(0.0;2);
    dstate[0] = -0.5*state[0];
    dstate[1] = -0.00001*state[1];
    
    dstate
}
```
