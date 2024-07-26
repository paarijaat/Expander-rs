use arith::MultiLinearPoly;
use crate::SumcheckMultilinearProdHelper;
use crate::GKRConfig;

#[derive(Clone, Debug, Default)]
pub struct GkrScratchpad<C: GKRConfig> {
    pub(crate) v_evals: Vec<C::Field>,
    pub(crate) hg_evals_5: Vec<C::Field>,
    pub(crate) hg_evals_1: Vec<C::Field>,

    pub(crate) eq_evals_at_rx: Vec<C::ChallengeField>,
    pub(crate) eq_evals_at_rz0: Vec<C::ChallengeField>,
    pub(crate) eq_evals_at_rz1: Vec<C::ChallengeField>,
    pub(crate) eq_evals_first_half: Vec<C::ChallengeField>,
    pub(crate) eq_evals_second_half: Vec<C::ChallengeField>,

    pub(crate) gate_exists_5: Vec<bool>,
    pub(crate) gate_exists_1: Vec<bool>,
}

impl<C: GKRConfig> GkrScratchpad<C> {
    pub(crate) fn new(max_num_input_var: usize, max_num_output_var: usize) -> Self {
        let max_input_num = 1 << max_num_input_var;
        let max_output_num = 1 << max_num_output_var;
        GkrScratchpad {
            v_evals: vec![C::Field::default(); max_input_num],
            hg_evals_5: vec![C::Field::default(); max_input_num],
            hg_evals_1: vec![C::Field::default(); max_input_num],

            eq_evals_at_rx: vec![C::ChallengeField::default(); max_input_num],
            eq_evals_at_rz0: vec![C::ChallengeField::default(); max_output_num],
            eq_evals_at_rz1: vec![C::ChallengeField::default(); max_output_num],
            eq_evals_first_half: vec![C::ChallengeField::default(); max_output_num],
            eq_evals_second_half: vec![C::ChallengeField::default(); max_output_num],

            gate_exists_5: vec![false; max_input_num],
            gate_exists_1: vec![false; max_input_num],
        }
    }
}


#[derive(Clone, Debug, Default)]
pub struct SumcheckMultilinearProdScratchpad<C: GKRConfig> {
    pub num_vars: usize,
    pub poly1: MultiLinearPoly<C::Field>,
    pub poly2: MultiLinearPoly<C::Field>,
    pub init_v: Vec<C::Field>,
    pub gate_exists: Vec<bool>,
    pub helper: SumcheckMultilinearProdHelper,
}

impl<C: GKRConfig> SumcheckMultilinearProdScratchpad<C> {
    pub fn new(poly1: &MultiLinearPoly<C::Field>, poly2: &MultiLinearPoly<C::Field>) -> Self {
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