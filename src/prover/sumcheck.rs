use arith::{Field, FieldSerde, VectorizedField};

use crate::{CircuitLayer, Config, GkrScratchpad, SumcheckMultilinearProdScratchpad, SumcheckGkrHelper, SumcheckMultilinearProdHelper, Transcript};
use arith::MultiLinearPoly;

// FIXME
#[allow(clippy::too_many_arguments)]
#[allow(clippy::type_complexity)]
pub fn sumcheck_prove_gkr_layer<F>(
    layer: &CircuitLayer<F>,
    rz0: &[Vec<F::BaseField>],
    rz1: &[Vec<F::BaseField>],
    alpha: &F::BaseField,
    beta: &F::BaseField,
    transcript: &mut Transcript,
    sp: &mut [GkrScratchpad<F>],
    config: &Config,
) -> (Vec<Vec<F::BaseField>>, Vec<Vec<F::BaseField>>)
where
    F: VectorizedField + FieldSerde,
    F::PackedBaseField: Field<BaseField = F::BaseField>,
{
    let mut helpers = vec![];
    assert_eq!(config.get_num_repetitions(), sp.len());
    for (j, sp_) in sp.iter_mut().enumerate() {
        helpers.push(SumcheckGkrHelper::new(
            layer, &rz0[j], &rz1[j], alpha, beta, sp_,
        ));
    }

    for i_var in 0..layer.input_var_num * 2 {
        for (j, helper) in helpers
            .iter_mut()
            .enumerate()
            .take(config.get_num_repetitions())
        {
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
            if j == 0 {
                log::trace!("i_var={} j={} evals: {:?} r: {:?}", i_var, j, evals, r);
            }
            helper.receive_challenge(i_var, r);
            if i_var == layer.input_var_num - 1 {
                log::trace!("vx claim: {:?}", helper.vx_claim());
                transcript.append_f(helper.vx_claim());
            }
        }
    }

    for (j, helper) in helpers
        .iter()
        .enumerate()
        .take(config.get_num_repetitions())
    {
        log::trace!("claimed vy[{}] = {:?}", j, helper.vy_claim());
        transcript.append_f(helper.vy_claim());
    }

    let rz0s = (0..config.get_num_repetitions())
        .map(|j| helpers[j].rx.clone()) // FIXME: clone might be avoided
        .collect();
    let rz1s = (0..config.get_num_repetitions())
        .map(|j| helpers[j].ry.clone()) // FIXME: clone might be avoided
        .collect();
    (rz0s, rz1s)
}



pub fn sumcheck_multilinear_prod<F>(
    transcript: &mut Transcript,
    sp: &mut SumcheckMultilinearProdScratchpad<F>,
) -> (Vec<F::BaseField>, (F,F))
where
    F: VectorizedField + FieldSerde,
    F::PackedBaseField: Field<BaseField = F::BaseField>,
{
    let mut randomness_sumcheck = Vec::<F::BaseField>::new();
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
    use arith::VectorizedFr;
    type F = VectorizedFr;

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

        let poly1 = MultiLinearPoly {
            var_num: num_vars,
            evals: evals1.clone()
        };

        let evals2 = vec![
            F::from(1 as u32),  // f(0,0)
            F::from(1 as u32),  // f(0,1)
            F::from(1 as u32), // f(1,0)
            F::from(1 as u32)  // f(1,1)
        ];

        let poly2 = MultiLinearPoly {
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