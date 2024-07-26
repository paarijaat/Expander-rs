use std::vec;

use log::info;
// use ark_std::{start_timer, end_timer};

use arith::MultiLinearPoly;
use expander_rs::{
    Transcript, sumcheck_multilinear_prod, SumcheckMultilinearProdScratchpad,
    BN254Config, M31ExtConfig, GKRConfig
};

#[macro_export]
macro_rules! log2i {
    ($i: expr) => {
        ($i as f32).log2().ceil() as usize
    };
}

fn sumcheck_multilinear_prod_test<C: GKRConfig>() {
    let num_vars: usize = 2;
    //let config = Config::bn254_config();

    let evals1 = vec![
        C::Field::from(4 as u32),  // f(0,0)
        C::Field::from(9 as u32),  // f(0,1)
        C::Field::from(16 as u32), // f(1,0)
        C::Field::from(25 as u32)  // f(1,1)
    ];

    let poly1 = MultiLinearPoly::<C::Field> {
        var_num: num_vars,
        evals: evals1.clone()
    };

    let evals2 = vec![
        C::Field::from(1 as u32),  // f(0,0)
        C::Field::from(1 as u32),  // f(0,1)
        C::Field::from(1 as u32), // f(1,0)
        C::Field::from(1 as u32)  // f(1,1)
    ];

    let poly2 = MultiLinearPoly::<C::Field> {
        var_num: num_vars,
        evals: evals2.clone()
    };

    let mut sp = SumcheckMultilinearProdScratchpad::<C>::new(&poly1, &poly2);
    let mut tp = Transcript::new();

    let (randomness_sumcheck, (p1,p2)) = sumcheck_multilinear_prod(
        &mut tp, 
        &mut sp,
    );

    let v1 = MultiLinearPoly::<C::Field>::eval_multilinear(
        &evals1,
        &randomness_sumcheck
    );

    let v2 = MultiLinearPoly::<C::Field>::eval_multilinear(
        &evals2,
        &randomness_sumcheck
    );

    info!("{:?}, {:?}", v1, v2);
    info!("{:?}, {:?}", p1, p2);
    assert_eq!(p1, v1);
    assert_eq!(p2, v2);
}

fn simple_tests<C: GKRConfig>() {
    let field_val = C::Field::from(5 as u32);
    info!("{:?}", field_val);
    
    let evals = vec![
        C::Field::from(4 as u32),  // f(0,0)
        C::Field::from(9 as u32),  // f(0,1)
        C::Field::from(16 as u32), // f(1,0)
        C::Field::from(25 as u32)  // f(1,1)
    ];
    //info!("evals: {:?}", evals);
    let eval_pt = vec![
        C::ChallengeField::from(11 as u32), // LSB
        C::ChallengeField::from(10 as u32)  // MSB
    ];

    let v = MultiLinearPoly::<C::Field>::eval_multilinear(
        &evals,
        eval_pt.as_slice()
    );

    info!("full eval {:?}", v);

    let eval_pt = vec![
        C::ChallengeField::from(0 as u32), // LSB
        C::ChallengeField::from(0 as u32)  // MSB
    ];

    let v = MultiLinearPoly::<C::Field>::eval_multilinear(
        &evals,
        eval_pt.as_slice()
    );

    info!("f(0,0) = {:?}", v);

    let v = MultiLinearPoly::<C::Field>::eval_multilinear(
        &evals,
        &[
            C::ChallengeField::from(1 as u32),  // LSB
            C::ChallengeField::from(0 as u32)   // MSB
        ],
    );

    info!("f(0,1) = {:?}", v);

    let v = MultiLinearPoly::<C::Field>::eval_multilinear(
        &evals,
        &[
            C::ChallengeField::from(0 as u32), 
            C::ChallengeField::from(1 as u32)
        ],
    );

    info!("f(1,0) = {:?}", v);

    let v = MultiLinearPoly::<C::Field>::eval_multilinear(
        &evals,
        &[
            C::ChallengeField::from(1 as u32), 
            C::ChallengeField::from(1 as u32)
        ],
    );

    info!("f(1,1) = {:?}", v);


    let poly = MultiLinearPoly {
        var_num: log2i!(evals.len()),
        evals: evals.clone()
    };

    // Fix LSB
    let poly_fix_var1 = poly.fix_variables_multilinear(&[C::ChallengeField::from(0 as u32)]);
    info!("poly_fix_var var_num = {:?}", &poly_fix_var1.var_num);
    info!("poly_fix_var evals = {:?}", &poly_fix_var1.evals);

    // Fix MSB
    let poly_fix_var2 = poly_fix_var1.fix_variables_multilinear(&[C::ChallengeField::from(1 as u32)]);
    info!("poly_fix_var var_num = {:?}", &poly_fix_var2.var_num);
    info!("poly_fix_var evals = {:?}", &poly_fix_var2.evals);

    // info!("evals: {:?}", evals);

}

fn main() {
    env_logger::init();

    simple_tests::<BN254Config>();
    sumcheck_multilinear_prod_test::<BN254Config>();

    // simple_tests::<M31ExtConfig>();
    // sumcheck_multilinear_prod_test::<M31ExtConfig>();
}
