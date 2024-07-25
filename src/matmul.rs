use std::vec;

use log::info;
// use ark_std::{start_timer, end_timer};

use arith::{Bn254DummyExt3, MultiLinearPoly, SimdM31Ext3, SimdField};
use expander_rs::{Transcript, sumcheck_multilinear_prod, SumcheckMultilinearProdScratchpad};

type F = Bn254DummyExt3;
//type F = SimdM31Ext3;


#[macro_export]
macro_rules! log2i {
    ($i: expr) => {
        ($i as f32).log2().ceil() as usize
    };
}

fn sumcheck_multilinear_prod_test() {
    let num_vars: usize = 2;
    //let config = Config::bn254_config();

    let evals1 = vec![
        F::from(4 as u32),  // f(0,0)
        F::from(9 as u32),  // f(0,1)
        F::from(16 as u32), // f(1,0)
        F::from(25 as u32)  // f(1,1)
    ];

    let poly1 = MultiLinearPoly::<F> {
        var_num: num_vars,
        evals: evals1.clone()
    };

    let evals2 = vec![
        F::from(1 as u32),  // f(0,0)
        F::from(1 as u32),  // f(0,1)
        F::from(1 as u32), // f(1,0)
        F::from(1 as u32)  // f(1,1)
    ];

    let poly2 = MultiLinearPoly::<F> {
        var_num: num_vars,
        evals: evals2.clone()
    };

    let mut sp = SumcheckMultilinearProdScratchpad::<F>::new(&poly1, &poly2);
    let mut tp = Transcript::new();

    let (randomness_sumcheck, (p1,p2)) = sumcheck_multilinear_prod(
        &mut tp, 
        &mut sp,
    );

    let v1 = MultiLinearPoly::<F>::eval_multilinear(
        &evals1,
        &randomness_sumcheck
    );

    let v2 = MultiLinearPoly::<F>::eval_multilinear(
        &evals2,
        &randomness_sumcheck
    );

    info!("{:?}, {:?}", v1, v2);
    info!("{:?}, {:?}", p1, p2);
    assert_eq!(p1, v1);
    assert_eq!(p2, v2);
}

fn main() {
    env_logger::init();

    let field_val = F::from(5 as u32);
    info!("{:?}", field_val);
    
    let evals = vec![
        F::from(4 as u32),  // f(0,0)
        F::from(9 as u32),  // f(0,1)
        F::from(16 as u32), // f(1,0)
        F::from(25 as u32)  // f(1,1)
    ];
    //info!("evals: {:?}", evals);
    let eval_pt = vec![
        <F as SimdField>::Scalar::from(11 as u32), // LSB
        <F as SimdField>::Scalar::from(10 as u32)  // MSB
    ];

    let v = MultiLinearPoly::<F>::eval_multilinear(
        &evals,
        eval_pt.as_slice()
    );

    info!("full eval {:?}", v);

    let eval_pt = vec![
        <F as SimdField>::Scalar::from(0 as u32), // LSB
        <F as SimdField>::Scalar::from(0 as u32)  // MSB
    ];

    let v = MultiLinearPoly::<F>::eval_multilinear(
        &evals,
        eval_pt.as_slice()
    );

    info!("f(0,0) = {:?}", v);

    let v = MultiLinearPoly::<F>::eval_multilinear(
        &evals,
        &[
            <F as SimdField>::Scalar::from(1 as u32),  // LSB
            <F as SimdField>::Scalar::from(0 as u32)   // MSB
        ],
    );

    info!("f(0,1) = {:?}", v);

    let v = MultiLinearPoly::<F>::eval_multilinear(
        &evals,
        &[
            <F as SimdField>::Scalar::from(0 as u32), 
            <F as SimdField>::Scalar::from(1 as u32)
        ],
    );

    info!("f(1,0) = {:?}", v);

    let v = MultiLinearPoly::<F>::eval_multilinear(
        &evals,
        &[
            <F as SimdField>::Scalar::from(1 as u32), 
            <F as SimdField>::Scalar::from(1 as u32)
        ],
    );

    info!("f(1,1) = {:?}", v);


    let poly = MultiLinearPoly {
        var_num: log2i!(evals.len()),
        evals: evals.clone()
    };

    // Fix LSB
    let poly_fix_var1 = poly.fix_variables_multilinear(&[<F as SimdField>::Scalar::from(0 as u32)]);
    info!("poly_fix_var var_num = {:?}", &poly_fix_var1.var_num);
    info!("poly_fix_var evals = {:?}", &poly_fix_var1.evals);

    // Fix MSB
    let poly_fix_var2 = poly_fix_var1.fix_variables_multilinear(&[<F as SimdField>::Scalar::from(1 as u32)]);
    info!("poly_fix_var var_num = {:?}", &poly_fix_var2.var_num);
    info!("poly_fix_var evals = {:?}", &poly_fix_var2.evals);

    // info!("evals: {:?}", evals);

    sumcheck_multilinear_prod_test();

}
