use std::{vec,iter::zip};
use log::info;
use ark_std::{start_timer, end_timer};

use arith::MultiLinearPoly;
use expander_rs::{
    Transcript, sumcheck_multilinear_prod, SumcheckMultilinearProdScratchpad,
    GKRConfig, BN254Config, Verifier,
};

#[macro_export]
macro_rules! log2i {
    ($i: expr) => {
        ($i as f32).log2().ceil() as usize
    };
}

fn sumcheck_multilinear_prod_test<C: GKRConfig>(num_vars: usize) {
    info!("========== Sumcheck {} vars =========", num_vars);
    //let config = Config::bn254_config();

    let mut evals1 = Vec::<C::Field>::new();
    for i in 0..(1 << num_vars) as u32 {
        evals1.push(C::Field::from(i))
    }
    // let evals1 = vec![
    //     C::Field::from(4 as u32),  // f(0,0)
    //     C::Field::from(9 as u32),  // f(0,1)
    //     C::Field::from(16 as u32), // f(1,0)
    //     C::Field::from(25 as u32)  // f(1,1)
    // ];

    let poly1 = MultiLinearPoly::<C::Field> {
        var_num: num_vars,
        evals: evals1.clone()
    };

    let mut evals2 = Vec::<C::Field>::new();
    for i in 0..(1 << num_vars) as u32 {
        evals2.push(C::Field::from(i+1))
    }
    // let evals2 = vec![
    //     C::Field::from(1 as u32),  // f(0,0)
    //     C::Field::from(1 as u32),  // f(0,1)
    //     C::Field::from(1 as u32), // f(1,0)
    //     C::Field::from(1 as u32)  // f(1,1)
    // ];

    let poly2 = MultiLinearPoly::<C::Field> {
        var_num: num_vars,
        evals: evals2.clone()
    };

    let mut sp = SumcheckMultilinearProdScratchpad::<C>::new(&poly1, &poly2);
    let mut tp = Transcript::new();

    let prove_time = start_timer!(|| format!("Sumcheck {} vars", num_vars));

    let (randomness_sumcheck, claimed_evals) = sumcheck_multilinear_prod(
        &mut tp, 
        &mut sp,
    );

    end_timer!(prove_time);

    let v1 = MultiLinearPoly::<C::Field>::eval_multilinear(
        &evals1,
        &randomness_sumcheck
    );

    let v2 = MultiLinearPoly::<C::Field>::eval_multilinear(
        &evals2,
        &randomness_sumcheck
    );

    let sum: C::Field = zip(&evals1, &evals2).map(|(a,b)| *a * *b).sum();
    info!("Sum:  {:?}", sum);
    //info!("Randomness:  {:?}", randomness_sumcheck);
    info!("Computed ev: {:?}, {:?}", v1, v2);
    info!("Claimed ev:  {:?}, {:?}", claimed_evals[0], claimed_evals[1]);
    assert_eq!(claimed_evals[0], v1);
    assert_eq!(claimed_evals[1], v2);

    let mut verified = false;
    let verifier = Verifier::<C>::default();
    verifier.verify_sumcheck(num_vars, &sum, &claimed_evals, &mut tp.proof, &mut verified);
    assert_eq!(verified, true);
    info!("Verified: true");

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

    info!("{}, {}, {}", 1 << 1, 1 << 2, 1 << 3);

    simple_tests::<BN254Config>();
    sumcheck_multilinear_prod_test::<BN254Config>(21);

    // simple_tests::<M31ExtConfig>();
    // sumcheck_multilinear_prod_test::<M31ExtConfig>();
}
