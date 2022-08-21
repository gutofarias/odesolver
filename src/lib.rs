use std::error::Error;
use std::fs::File;
use std::io::Write;

#[derive(Clone)]
pub struct EdoParam {
    pub time : f64,
    pub tend : f64,
    pub step : f64,
    pub relStepOut: usize,
}

pub mod solver_trait;
pub mod solver_vector;
