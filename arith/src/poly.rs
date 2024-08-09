use ark_std::{end_timer, start_timer};

use crate::{Field, SimdField};

#[derive(Debug, Clone, Default, PartialEq, Eq)]
/// Definition for an MLE, with an associated type F.
pub struct MultiLinearPoly<F: Field + SimdField> {
    /// Number of variables in an MLE
    pub var_num: usize,
    /// MLE Evaluations
    pub evals: Vec<F>,
}

impl<F: Field + SimdField> MultiLinearPoly<F> {
    pub fn eval_multilinear(evals: &[F], x: &[F::Scalar]) -> F {
        let timer = start_timer!(|| format!("eval mle with {} vars", x.len()));
        assert_eq!(1 << x.len(), evals.len());
        let mut scratch = evals.to_vec();
        let mut cur_eval_size = evals.len() >> 1;
        for r in x.iter() {
            log::trace!("scratch: {:?}", scratch);
            for i in 0..cur_eval_size {
                scratch[i] = scratch[i * 2] + (scratch[i * 2 + 1] - scratch[i * 2]).scale(r);
            }
            cur_eval_size >>= 1;
        }
        end_timer!(timer);
        scratch[0]
    }

    pub fn fix_variables_multilinear_lsb_first(&self, partial_point: &[F::Scalar]) -> Self {
        let timer = start_timer!(|| format!("fix variable mle with {} vars", partial_point.len()));
        assert!(partial_point.len() <= self.var_num, "invalid size of partial point");
        let mut scratch = self.evals.to_vec();
        log::trace!("scratch({}): {:?}", scratch.len(), scratch[..scratch.len()].to_vec());
        let mut cur_eval_size = scratch.len() >> 1;
        for r in partial_point.iter() {
            log::trace!("Fix lsb to = {:?}", r);
            for i in 0..cur_eval_size {
                scratch[i] = scratch[i * 2] + (scratch[i * 2 + 1] - scratch[i * 2]).scale(r);
            }
            log::trace!("scratch({}): {:?}", cur_eval_size, scratch[..cur_eval_size].to_vec());
            cur_eval_size >>= 1;
        }
        let num_remaining_vars = self.var_num - partial_point.len();
        end_timer!(timer);
        Self {
            var_num: num_remaining_vars,
            evals: scratch[..(1 << num_remaining_vars)].to_vec()
        }
    }

    pub fn fix_variables_multilinear_msb_first(&self, partial_point: &[F::Scalar]) -> Self {
        let timer = start_timer!(|| format!("fix variable mle with {} vars", partial_point.len()));
        assert!(partial_point.len() <= self.var_num, "invalid size of partial point");
        let mut scratch = self.evals.to_vec();
        log::trace!("scratch({}): {:?}", scratch.len(), scratch[..scratch.len()].to_vec());
        let mut cur_eval_size = scratch.len() >> 1;
        for r in partial_point.iter() {
            log::trace!("Fix msb to = {:?}", r);
            for i in 0..cur_eval_size {
                scratch[i] = scratch[i] + (scratch[i + cur_eval_size] - scratch[i]).scale(r);
            }
            log::trace!("scratch({}): {:?}", cur_eval_size, scratch[..cur_eval_size].to_vec());
            cur_eval_size >>= 1;
        }
        let num_remaining_vars = self.var_num - partial_point.len();
        end_timer!(timer);
        Self {
            var_num: num_remaining_vars,
            evals: scratch[..(1 << num_remaining_vars)].to_vec()
        }
    }
}

