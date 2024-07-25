use arith::{Field, SimdField, MultiLinearPoly};
use crate::SumcheckMultilinearProdHelper;

#[derive(Clone, Debug, Default)]
pub struct GkrScratchpad<F: Field + SimdField> {
    pub(crate) v_evals: Vec<F>,
    pub(crate) hg_evals: Vec<F>,

    pub(crate) eq_evals_at_rx: Vec<F::Scalar>,
    pub(crate) eq_evals_at_rz0: Vec<F::Scalar>,
    pub(crate) eq_evals_at_rz1: Vec<F::Scalar>,
    pub(crate) eq_evals_first_half: Vec<F::Scalar>,
    pub(crate) eq_evals_second_half: Vec<F::Scalar>,

    pub(crate) gate_exists: Vec<bool>,
}

impl<F: Field + SimdField> GkrScratchpad<F> {
    pub(crate) fn new(max_num_input_var: usize, max_num_output_var: usize) -> Self {
        let max_input_num = 1 << max_num_input_var;
        let max_output_num = 1 << max_num_output_var;
        GkrScratchpad {
            v_evals: vec![F::default(); max_input_num],
            hg_evals: vec![F::default(); max_input_num],

            eq_evals_at_rx: vec![F::Scalar::default(); max_input_num],
            eq_evals_at_rz0: vec![F::Scalar::default(); max_output_num],
            eq_evals_at_rz1: vec![F::Scalar::default(); max_output_num],
            eq_evals_first_half: vec![F::Scalar::default(); max_output_num],
            eq_evals_second_half: vec![F::Scalar::default(); max_output_num],

            gate_exists: vec![false; max_input_num],
        }
    }
}


#[derive(Clone, Debug, Default)]
pub struct SumcheckMultilinearProdScratchpad<F: Field + SimdField> {
    pub num_vars: usize,
    pub poly1: MultiLinearPoly<F>,
    pub poly2: MultiLinearPoly<F>,
    pub init_v: Vec<F>,
    pub gate_exists: Vec<bool>,
    pub helper: SumcheckMultilinearProdHelper,
}

impl<F: Field + SimdField> SumcheckMultilinearProdScratchpad<F> {
    pub fn new(poly1: &MultiLinearPoly<F>, poly2: &MultiLinearPoly<F>) -> Self {
        let num_vars = poly1.var_num;
        let num_evals = 1 << num_vars;
        assert_eq!(num_evals, poly1.evals.len());
        assert_eq!(num_evals, poly2.evals.len());
        SumcheckMultilinearProdScratchpad {
            num_vars,
            poly1: poly1.clone(),
            poly2: poly2.clone(),
            init_v: poly1.evals.clone(),
            gate_exists: vec![true; num_evals],
            helper: SumcheckMultilinearProdHelper::new(num_vars),
        }
    }
}