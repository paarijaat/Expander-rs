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

    pub fn fix_variables_multilinear(&self, partial_point: &[F::Scalar]) -> Self {
        let timer = start_timer!(|| format!("fix variable mle with {} vars", partial_point.len()));
        assert!(partial_point.len() <= self.var_num, "invalid size of partial point");
        let mut scratch = self.evals.to_vec();
        let mut cur_eval_size = scratch.len() >> 1;
        for r in partial_point.iter() {
            log::trace!("scratch: {:?}", scratch);
            for i in 0..cur_eval_size {
                scratch[i] = scratch[i * 2] + (scratch[i * 2 + 1] - scratch[i * 2]).scale(r);
            }
            cur_eval_size >>= 1;
        }
        let num_remaining_vars = self.var_num - partial_point.len();
        end_timer!(timer);
        Self {
            var_num: num_remaining_vars,
            evals: scratch[..(1 << num_remaining_vars)].to_vec()
        }
    }

    // pub fn fix_variables(&self, partial_point: &[F::Scalar]) -> Self {
    //     assert!(
    //         partial_point.len() <= self.var_num,
    //         "invalid size of partial point"
    //     );

    //     let mut poly = self.evals.to_vec();
    //     let nv = self.var_num;
    //     let dim = partial_point.len();
    //     // evaluate single variable of partial point from left to right
    //     for i in 1..dim + 1 {
    //         let r = partial_point[i - 1];
    //         for b in 0..(1 << (nv - i)) {
    //             let left = poly[b << 1];
    //             let right = poly[(b << 1) + 1];
    //             poly[b] = left + (right - left).scale(&r);
    //         }
    //     }
    //     Self {
    //         var_num: nv - dim,
    //         evals: poly[..(1 << (nv - dim))].to_vec()
    //     }
    // }
}

