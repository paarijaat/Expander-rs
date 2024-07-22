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
    m1: &mut MultiLinearPoly<F>,
    m2: &mut MultiLinearPoly<F>,
    transcript: &mut Transcript,
    sp: &mut SumcheckMultilinearProdScratchpad<F>,
    init_v: &[F],
) -> (Vec<F::BaseField>, (F,F))
where
    F: VectorizedField + FieldSerde,
    F::PackedBaseField: Field<BaseField = F::BaseField>,
{
    let mut randomness_during_sumcheck = Vec::<F::BaseField>::new();
    let mut claimed_evals_m1_m2 = (F::zero(), F::zero());
    let mut helper = SumcheckMultilinearProdHelper::new(m1.var_num);

    for i_var in 0..m1.var_num {
        /* 
        fn poly_eval_at<F: VectorizedField>(
            &self,
            var_idx: usize,
            degree: usize,
            bk_f: &mut [F],
            bk_hg: &mut [F],
            init_v: &[F],
            gate_exists: &[bool],
        ) -> [F; 3]
        */
        let evals = helper.poly_eval_at::<F>(
            i_var, 
            2, 
            m1.evals.as_mut_slice(),
            m2.evals.as_mut_slice(),
            init_v, // TODO fix this
            &sp.gate_exists
        );

        transcript.append_f(evals[0]);
        transcript.append_f(evals[1]);
        transcript.append_f(evals[2]);

        let r = transcript.challenge_f::<F>();
        randomness_during_sumcheck.push(r.clone());
        /*
            fn receive_challenge<F: VectorizedField>(
                &mut self,
                var_idx: usize,
                r: F::BaseField,
                bk_f: &mut [F],
                bk_hg: &mut [F],
                init_v: &[F],
                gate_exists: &mut [bool],
            ) where
                F::PackedBaseField: Field<BaseField = F::BaseField>,
        */
        helper.receive_challenge::<F>(
            i_var, 
            r, 
            m1.evals.as_mut_slice(), 
            m2.evals.as_mut_slice(),
            init_v, // TODO fix this
            &mut sp.gate_exists
        );
    }

    // Claimed evaluations of m1 and m2
    // log::trace!("vx claim: {:?}", helper.vx_claim());
    claimed_evals_m1_m2.0 = m1.evals[0].clone();
    transcript.append_f(m1.evals[0]);

    // log::trace!("claimed vy[{}] = {:?}", j, helper.vy_claim());
    claimed_evals_m1_m2.1 = m2.evals[0].clone();
    transcript.append_f(m2.evals[0]);

    (randomness_during_sumcheck, claimed_evals_m1_m2)
}


#[cfg(test)]
mod tests {
    use super::*;
    use arith::{Field, VectorizedFr};
    type F = VectorizedFr;

    #[test]
    fn sumcheck_multilinear_prod_test() {
        let num_vars: usize = 2;
        //let config = Config::bn254_config();

        let mut sp = SumcheckMultilinearProdScratchpad::<F>::new(num_vars);
        let mut tp = Transcript::new();

        let evals1 = vec![
            F::from(4 as u32),  // f(0,0)
            F::from(9 as u32),  // f(0,1)
            F::from(16 as u32), // f(1,0)
            F::from(25 as u32)  // f(1,1)
        ];

        let mut poly1 = MultiLinearPoly {
            var_num: num_vars,
            evals: evals1.clone()
        };

        let evals2 = vec![
            F::from(1 as u32),  // f(0,0)
            F::from(1 as u32),  // f(0,1)
            F::from(1 as u32), // f(1,0)
            F::from(1 as u32)  // f(1,1)
        ];

        let mut poly2 = MultiLinearPoly {
            var_num: num_vars,
            evals: evals2.clone()
        };

        // let mut p1 = F::zero();
        // let mut p2 = F::zero();

        let (randomness_sumcheck, (p1,p2)) = sumcheck_multilinear_prod(
            &mut poly1, 
            &mut poly2,
            &mut tp, 
            &mut sp,
            &evals1.clone()
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