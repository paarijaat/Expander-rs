use arith::{SimdField, Field, FieldSerde};
use crate::{
    CircuitLayer, Config, GkrScratchpad, SumcheckGkrHelper, SumcheckGkrSquareHelper, Transcript, SumcheckMultilinearProdScratchpad
};

// FIXME
#[allow(clippy::too_many_arguments)]
#[allow(clippy::type_complexity)]
pub fn sumcheck_prove_gkr_layer<F>(
    layer: &CircuitLayer<F>,
    rz0: &[F::Scalar],
    rz1: &[F::Scalar],
    alpha: &F::Scalar,
    beta: &F::Scalar,
    transcript: &mut Transcript,
    sp: &mut GkrScratchpad<F>,
    _config: &Config,
) -> (Vec<F::Scalar>, Vec<F::Scalar>)
where
    F: Field + FieldSerde + SimdField,
{
    let mut helper = SumcheckGkrHelper::new(layer, rz0, rz1, alpha, beta, sp);

    for i_var in 0..layer.input_var_num * 2 {
        if i_var == 0 {
            helper.prepare_g_x_vals()
        }

        if i_var == layer.input_var_num {
            let vx_claim = helper.vx_claim();
            helper.prepare_h_y_vals(vx_claim)
        }

        let evals = helper.poly_evals_at(i_var, 2);

        transcript.append_f(evals[0]);
        transcript.append_f(evals[1]);
        transcript.append_f(evals[2]);

        let r = transcript.challenge_f::<F>();

        log::trace!("i_var={} evals: {:?} r: {:?}", i_var, evals, r);

        helper.receive_challenge(i_var, r);
        if i_var == layer.input_var_num - 1 {
            log::trace!("vx claim: {:?}", helper.vx_claim());
            transcript.append_f(helper.vx_claim());
        }
    }

    log::trace!("claimed vy = {:?}", helper.vy_claim());
    transcript.append_f(helper.vy_claim());

    let rz0 = helper.rx.clone();
    let rz1 = helper.ry.clone();
    (rz0, rz1)
}

// FIXME
#[allow(clippy::too_many_arguments)]
#[allow(clippy::type_complexity)]
#[allow(clippy::needless_range_loop)]
pub fn sumcheck_prove_gkr_square_layer<F>(
    layer: &CircuitLayer<F>,
    rz0: &[F::Scalar],
    transcript: &mut Transcript,
    sp: &mut GkrScratchpad<F>,
    _config: &Config,
) -> Vec<F::Scalar>
where
    F: Field + FieldSerde + SimdField,
{
    const D: usize = 7;
    let mut helper = SumcheckGkrSquareHelper::new(layer, rz0, sp);

    for i_var in 0..layer.input_var_num {
        if i_var == 0 {
            helper.prepare_g_x_vals()
        }

        let evals: [F; D] = helper.poly_evals_at(i_var);

        for deg in 0..D {
            transcript.append_f(evals[deg]);
        }

        let r = transcript.challenge_f::<F>();

        log::trace!("i_var={} evals: {:?} r: {:?}", i_var, evals, r);

        helper.receive_challenge(i_var, r);
    }

    log::trace!("claimed vx = {:?}", helper.vx_claim());
    transcript.append_f(helper.vx_claim());

    helper.rx.clone()
}



pub fn sumcheck_multilinear_prod<F>(
    transcript: &mut Transcript,
    sp: &mut SumcheckMultilinearProdScratchpad<F>,
) -> (Vec<F::Scalar>, (F,F))
where
    F: Field + FieldSerde + SimdField,
{
    let mut randomness_sumcheck = Vec::<F::Scalar>::new();
    let mut claimed_evals_m1_m2 = (F::zero(), F::zero());

    for i_var in 0..sp.num_vars {
        // Computes the three values (evaluations of the poly for the verifier)
        let evals = sp.helper.poly_eval_at::<F>(
            i_var, 
            2, 
            sp.poly1.evals.as_mut_slice(),
            sp.poly2.evals.as_mut_slice(),
            &sp.init_v, // TODO: fix this, why do we need to send this
            &sp.gate_exists  // TODO: fix this, why do we need to send this
        );

        // Append the poly sent to verifier to the transcript
        transcript.append_f(evals[0]);
        transcript.append_f(evals[1]);
        transcript.append_f(evals[2]);

        // Create the next randomness (fiat-shamir)
        let r = transcript.challenge_f::<F>();
        randomness_sumcheck.push(r.clone());

        // Fix the next variable using the fiat-shamir randomness
        sp.helper.receive_challenge::<F>(
            i_var, 
            r, 
            sp.poly1.evals.as_mut_slice(), 
            sp.poly2.evals.as_mut_slice(),
            &sp.init_v, // TODO: fix this, why do we need to send this
            &mut sp.gate_exists  // TODO: fix this, why do we need to send this
        );
    }

    // Claimed evaluations of m1 and m2
    // log::trace!("vx claim: {:?}", helper.vx_claim());
    claimed_evals_m1_m2.0 = sp.poly1.evals[0].clone();
    transcript.append_f(claimed_evals_m1_m2.0);

    // log::trace!("claimed vy[{}] = {:?}", j, helper.vy_claim());
    claimed_evals_m1_m2.1 = sp.poly2.evals[0].clone();
    transcript.append_f(claimed_evals_m1_m2.1);

    (randomness_sumcheck, claimed_evals_m1_m2)
}


#[cfg(test)]
mod tests {
    use super::*;
    use arith::{MultiLinearPoly, Bn254DummyExt3};
    type F = Bn254DummyExt3;

    #[test]
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

        println!("{:?}, {:?}", v1, v2);
        println!("{:?}, {:?}", p1, p2);
        assert_eq!(p1, v1);
        assert_eq!(p2, v2);
    }
}