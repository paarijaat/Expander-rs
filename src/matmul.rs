use std::vec;

use halo2curves::bn256::Fr;

use log::info;
// use ark_std::{start_timer, end_timer};
use arith::MultiLinearPoly;


pub mod circuit;
pub use circuit::*;

pub mod config;
pub use config::*;

pub mod hash;
pub use hash::*;

pub mod poly_commit;
pub use poly_commit::*;

pub mod prover;
pub use prover::*;

// type F = VectorizedFr;
// type Fbase = <VectorizedFr as Field>::BaseField;
type F = Fr;
// type Fbase = <VectorizedFr as Field>::BaseField;



#[macro_export]
macro_rules! log2i {
    ($i: expr) => {
        ($i as f32).log2().ceil() as usize
    };
}

fn main() {
    env_logger::init();
    
    // input_vals: MultiLinearPoly::<F> {
    //     var_num: layer_seg.i_var_num,
    //     evals: vec![],
    // },

    let field_val = F::from(5 as u32);
    let base_field_val = F::from(10 as u32);
    // info!("{:?}", field_val);
    // info!("{:?}", base_field_val);
    // info!("");
    
    let evals = vec![
        F::from(4 as u32),  // f(0,0)
        F::from(9 as u32),  // f(0,1)
        F::from(16 as u32), // f(1,0)
        F::from(25 as u32)  // f(1,1)
    ];
    //info!("evals: {:?}", evals);
    let v = MultiLinearPoly::<F>::eval_multilinear(
        &evals,
        &[
            F::from(11 as u32), // LSB
            F::from(10 as u32)  // MSB
        ],
    );

    info!("full eval {:?}", v);

    let v = MultiLinearPoly::<F>::eval_multilinear(
        &evals,
        &[
            F::from(0 as u32), 
            F::from(0 as u32)
        ],
    );

    info!("f(0,0) = {:?}", v);

    let v = MultiLinearPoly::<F>::eval_multilinear(
        &evals,
        &[
            F::from(1 as u32),  // LSB
            F::from(0 as u32)   // MSB
        ],
    );

    info!("f(0,1) = {:?}", v);

    let v = MultiLinearPoly::<F>::eval_multilinear(
        &evals,
        &[
            F::from(0 as u32), 
            F::from(1 as u32)
        ],
    );

    info!("f(1,0) = {:?}", v);

    let v = MultiLinearPoly::<F>::eval_multilinear(
        &evals,
        &[
            F::from(1 as u32), 
            F::from(1 as u32)
        ],
    );

    info!("f(1,1) = {:?}", v);


    let poly = MultiLinearPoly {
        var_num: log2i!(evals.len()),
        evals: evals.clone()
    };

    // Fix LSB
    let poly_fix_var1 = poly.fix_variables(&[F::from(0 as u32)]);
    info!("poly_fix_var var_num = {:?}", &poly_fix_var1.var_num);
    info!("poly_fix_var evals = {:?}", &poly_fix_var1.evals);

    // Fix MSB
    let poly_fix_var2 = poly_fix_var1.fix_variables(&[F::from(1 as u32)]);
    info!("poly_fix_var var_num = {:?}", &poly_fix_var2.var_num);
    info!("poly_fix_var evals = {:?}", &poly_fix_var2.evals);
    // eq_eval_at()

    // eq_evals_at_primitive(
    //     &[F::from(2 as u32)], 
    //     &F::one(), 
    //     &mut evals
    // );

    // info!("evals: {:?}", evals);


    // eq_evals_at_primitive(
    //     &[F::from(3 as u32)], 
    //     &F::one(), 
    //     &mut evals
    // );

    // info!("evals: {:?}", evals);

}
