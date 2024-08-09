use std::{iter::zip, vec};
use arith::Field;
use log::{info, debug};
use ark_std::{end_timer, start_timer};

use arith::MultiLinearPoly;
use expander_rs::{
    Transcript, sumcheck_multilinear_prod, SumcheckMultilinearProdScratchpad,
    GKRConfig, BN254Config, Verifier,
};
use csv;
use csv::StringRecord;

#[macro_export]
macro_rules! log2i {
    ($i: expr) => {
        ($i as f32).log2().ceil() as usize
    };
}

// Computes the binary reverse number of n
// with zeros padding to the left to N bits i.e.
// if n = 7 and N = 4 then bin(n) = 111, padded(bin(n),N) = 0111
// and reversed(padded(bin(n))) = 1110 = 14
// Used to determine the index of matrix values inside the MLE
#[macro_export]
macro_rules! rbin {
    ($n: expr, $nb_bits: expr) => {
        if $nb_bits == 0{
            0
        } else {
        ($n).reverse_bits() >> ((usize::BITS as usize) - $nb_bits)
        }
    };
}

pub type Matrix<T> = Vec<Vec<T>>;

fn read_matrix<C: GKRConfig>(filename: &str) -> Matrix<C::Field> {

    let mut reader = csv::ReaderBuilder::new()
    .has_headers(false)
    .from_path(filename).unwrap();
    // let mut reader = csv::Reader::from_path(filename)?;
    let n_r = reader.records().count();

    let mut reader2 = csv::ReaderBuilder::new()
    .has_headers(false)
    .from_path(filename).unwrap();

    let start = reader2.position().clone();

    let mut record = StringRecord::new();

    reader2.read_record(&mut record).unwrap();
    let n_c = record.len();
        
    // println!("Rows:{}", n_r);
    // println!("Cols:{}", n_c);

    let mut mat = vec![vec![C::Field::zero(); n_c]; n_r];
    let mut cr = 0;

    reader2.seek(start).unwrap();
    for result in reader2.records() {
        let record = result.unwrap();
        for i in 0..record.len() {
            mat[cr][i] = C::Field::from((record.get(i).unwrap()).parse::<u32>().unwrap());
        }
        cr = cr + 1;
    }
    // println!("Check:{}", cr);
    return mat;
}


pub fn matrix_to_mle<C: GKRConfig>(
    matrix: &Matrix<C::Field>,
    transposed: bool,
    reverse_binary_index: bool,
) -> (MultiLinearPoly<C::Field>, usize, usize) {
    /*
    Take an input matrix and evaluate it as an MLE.
    Rows are concatenated in the MLE. Namely, the ordering of the elements is:
    [(0,0), (0,1), ..., (0, n_cols - 1), (1,0), ..., (n_rows - 1, n_cols - 1)]
    For n_cols and n_rows being a power of 2. Otherwise, we pad with zeroes.
    */
    let num_rows = matrix.len();
    let num_cols = matrix[0].len();
    let num_vars_rows = log2i!(num_rows);
    let num_vars_cols = log2i!(num_cols);

    let n = num_vars_rows + num_vars_cols;

    // Initialize an MLE of zeros (field elements)
    let mut poly_evals = vec![C::Field::zero() ; 1 << n];
    if transposed {
        for j in 0..num_cols {
            let index = j * (1 << num_vars_rows);
            for i in 0..num_rows {
                let mut index_to_use = index + i;
                if reverse_binary_index {
                    index_to_use = rbin!(index_to_use, n);
                    // let rev_index = rbin!(index + i, n);
                }
                poly_evals[index_to_use] = matrix[i][j];
            }
        }
    } else {
        for i in 0..num_rows {
            let index = i * (1 << num_vars_cols);
            for j in 0..num_cols {
                let mut index_to_use = index + j;
                if reverse_binary_index {
                    index_to_use = rbin!(index_to_use, n);
                    // let rev_index = rbin!(index + j, n);
                }
                poly_evals[index_to_use] = matrix[i][j];
            }
        }
    }

    let poly = MultiLinearPoly::<C::Field> {
        var_num: n,
        evals: poly_evals
    };
    return (poly, num_vars_rows, num_vars_cols);
}

fn get_dummy_mat<C: GKRConfig>(num_rows: usize, num_cols: usize) -> Matrix<C::Field> {
    let mut count: u32 = 1;
    let mut mat: Matrix<C::Field> = vec![vec![C::Field::zero(); num_cols]; num_rows];
    for row in 0..num_rows {
        for col in 0..num_cols {
            mat[row][col] = C::Field::from(count);
            count += 1;
        }
    }
    return mat
}

#[allow(dead_code)]
fn matrix_tests<C: GKRConfig>() {
    let num_rows: usize = 4;
    let num_cols: usize = 4;
    let mat = get_dummy_mat::<C>(num_rows, num_cols);

    let (mat1_mle, num_row_vars, num_col_vars) = matrix_to_mle::<C>(&mat, false, true);

    let num_vars: usize = num_row_vars + num_col_vars;

    let mut mat2_evals = Vec::<C::Field>::new();
    for i in 1..=(1 << num_vars) as u32 {
        mat2_evals.push(C::Field::from(i))
    }

    //info!("Mat: {:?}", mat_evals);
    let mat2_mle = MultiLinearPoly::<C::Field> {
        var_num: num_vars,
        evals: mat2_evals.clone()
    };

    let eval_pt: Vec<u32> = vec![0,1,0,0];

    let eval_pt_field =         &vec![
        C::ChallengeField::from(eval_pt[0]), 
        C::ChallengeField::from(eval_pt[1]),
        C::ChallengeField::from(eval_pt[2]),
        C::ChallengeField::from(eval_pt[3]),
    ];

    // let v1 = MultiLinearPoly::<C::Field>::eval_multilinear(
    //     &mat1_mle.evals,
    //     &eval_pt_field.as_slice()
    // );
    // info!("mat1_bin_reverse({:?}) = {:?}", eval_pt, v1);

    // let v2 = MultiLinearPoly::<C::Field>::eval_multilinear(
    //     &mat2_mle.evals,
    //     &eval_pt_field.as_slice()
    // );
    // info!("mat2_simple({:?}) = {:?}", eval_pt, v2);

    debug!("mat1_mle.fix_variables_multilinear_lsb_first");
    let mat1_mle_lsb_first_eval = mat1_mle.fix_variables_multilinear_lsb_first(eval_pt_field);
    debug!("mat1_mle_lsb_first_eval: {:?}", mat1_mle_lsb_first_eval.evals[0]);
    // info!("mat1_mle.fix_variables_multilinear_msb_first");
    // let _ = &mat1_mle.fix_variables_multilinear_msb_first(eval_pt_field);
    // info!("mat2_mle.fix_variables_multilinear_lsb_first");
    // let _ = &mat2_mle.fix_variables_multilinear_lsb_first(eval_pt_field);
    debug!("mat2_mle.fix_variables_multilinear_msb_first");
    let mat2_mle_msb_first_eval = mat2_mle.fix_variables_multilinear_msb_first(eval_pt_field);
    debug!("mat2_mle_msb_first_eval: {:?}", mat2_mle_msb_first_eval.evals[0]);

    let row_eval_pt =         &vec![
        C::ChallengeField::from(1), 
        C::ChallengeField::from(0),
    ];
    let picked_row_mle = mat1_mle.fix_variables_multilinear_lsb_first(row_eval_pt);
    debug!("picked_row_mle: {:?}", picked_row_mle.evals);

    let col_eval_pt =         &vec![
        C::ChallengeField::from(1), 
        C::ChallengeField::from(0),
    ];
    let picked_col_mle = mat2_mle.fix_variables_multilinear_lsb_first(col_eval_pt);
    debug!("picked_col_mle: {:?}", picked_col_mle.evals);

    debug!("=======================================");
    let (mat_mle_rev, _, _) = matrix_to_mle::<C>(&mat, false, true);
    debug!("mat_mle_rev: {:?}", mat_mle_rev.evals);
    let (mat_mle_norm, _, _) = matrix_to_mle::<C>(&mat, false, false);
    debug!("mat_mle_norm: {:?}", mat_mle_norm.evals);

    let fix_var_pt =         &vec![
        C::ChallengeField::from(1), 
        C::ChallengeField::from(0),
    ];

    let mat_mle_rev_picked = mat_mle_rev.fix_variables_multilinear_lsb_first(fix_var_pt);
    debug!("mat_mle_rev_picked: {:?}", mat_mle_rev_picked.evals);

    let mat_mle_norm_picked = mat_mle_norm.fix_variables_multilinear_lsb_first(fix_var_pt);
    debug!("mat_mle_norm picked: {:?}", mat_mle_norm_picked.evals);


    debug!("==============Not using reverse and no transpose==================");
    let (mat_mle_msb, _, _) = matrix_to_mle::<C>(&mat, false, false);
    debug!("mat_mle_msb: {:?}", mat_mle_msb.evals);
    let (mat_mle_lsb, _, _) = matrix_to_mle::<C>(&mat, false, false);
    debug!("mat_mle_lsb: {:?}", mat_mle_lsb.evals);

    let fix_var_pt =         &vec![
        C::ChallengeField::from(1), 
        C::ChallengeField::from(0),
    ];

    let mat_mle_msb_picked = mat_mle_msb.fix_variables_multilinear_msb_first(fix_var_pt);
    debug!("mat_mle_msb_picked: {:?}", mat_mle_msb_picked.evals);

    let mat_mle_lsb_picked = mat_mle_lsb.fix_variables_multilinear_lsb_first(fix_var_pt);
    debug!("mat_mle_lsb_picked: {:?}", mat_mle_lsb_picked.evals);

    let mut sp = SumcheckMultilinearProdScratchpad::<C>::new(&mat_mle_msb_picked, &mat_mle_lsb_picked);
    let mut tp = Transcript::new();

    let (randomness_sumcheck, claimed_evals) = sumcheck_multilinear_prod(
        &mut tp, 
        &mut sp,
    );

    let v1 = MultiLinearPoly::<C::Field>::eval_multilinear(
        &mat_mle_msb_picked.evals,
        &randomness_sumcheck
    );

    let v2 = MultiLinearPoly::<C::Field>::eval_multilinear(
        &mat_mle_lsb_picked.evals,
        &randomness_sumcheck
    );

    let sum: C::Field = zip(&mat_mle_msb_picked.evals, &mat_mle_lsb_picked.evals).map(|(a,b)| *a * *b).sum();
    info!("Sum:  {:?}", sum);
    //info!("Randomness:  {:?}", randomness_sumcheck);
    info!("Computed ev: {:?}, {:?}", v1, v2);
    info!("Claimed ev:  {:?}, {:?}", claimed_evals[0], claimed_evals[1]);
    assert_eq!(claimed_evals[0], v1);
    assert_eq!(claimed_evals[1], v2);

    let mut verified = false;
    let verifier = Verifier::<C>::default();
    verifier.verify_sumcheck(mat_mle_msb_picked.var_num, &sum, &claimed_evals, &mut tp.proof, &mut verified);
    assert_eq!(verified, true);
    info!("Verified: true");


    let mut randomness_for_msb_initiated_poly = fix_var_pt.clone();
    for r in randomness_sumcheck.iter().rev() {
        randomness_for_msb_initiated_poly.push(*r);
    }

    let mat_mle_msb_full_msb_randomness = mat_mle_msb.fix_variables_multilinear_msb_first(&randomness_for_msb_initiated_poly);
    debug!("mat_mle_msb_full_msb_randomness: {:?}", mat_mle_msb_full_msb_randomness.evals[0]);

    randomness_for_msb_initiated_poly.reverse();
    debug!("randomness_for_msb_initiated_poly: {:?}", randomness_for_msb_initiated_poly);
    
    let v1_from_orig_msb_poly = MultiLinearPoly::<C::Field>::eval_multilinear(
        &mat_mle_msb.evals,
        &randomness_for_msb_initiated_poly
    );
    assert_eq!(claimed_evals[0], v1_from_orig_msb_poly);
    assert_eq!(claimed_evals[0], mat_mle_msb_full_msb_randomness.evals[0]);


    let mut randomness_for_lsb_initiated_poly = fix_var_pt.clone();
    for r in randomness_sumcheck.iter() {
        randomness_for_lsb_initiated_poly.push(*r);
    }
    let v2_from_orig_lsb_poly = MultiLinearPoly::<C::Field>::eval_multilinear(
        &mat_mle_lsb.evals,
        &randomness_for_lsb_initiated_poly
    );
    debug!("randomness_for_lsb_initiated_poly: {:?}", randomness_for_lsb_initiated_poly);
    assert_eq!(claimed_evals[1], v2_from_orig_lsb_poly);
}


#[allow(dead_code)]
fn sumcheck_multilinear_prod_test<C: GKRConfig>(num_vars: usize) {
    info!("========== Sumcheck {} vars =========", num_vars);

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


#[allow(dead_code)]
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
    let poly_fix_var1 = poly.fix_variables_multilinear_lsb_first(&[C::ChallengeField::from(0 as u32)]);
    info!("poly_fix_var var_num = {:?}", &poly_fix_var1.var_num);
    info!("poly_fix_var evals = {:?}", &poly_fix_var1.evals);

    // Fix MSB
    let poly_fix_var2 = poly_fix_var1.fix_variables_multilinear_lsb_first(&[C::ChallengeField::from(1 as u32)]);
    info!("poly_fix_var var_num = {:?}", &poly_fix_var2.var_num);
    info!("poly_fix_var evals = {:?}", &poly_fix_var2.evals);

    // info!("evals: {:?}", evals);

}

fn main() {
    env_logger::init();

    info!("{}, {}, {}, {}, {}", 1 << 1, 1 << 2, 1 << 3, 5 >> 1, 5 >> 2);

    //simple_tests::<BN254Config>();
    // sumcheck_multilinear_prod_test::<BN254Config>(21);

    matrix_tests::<BN254Config>();

    // simple_tests::<M31ExtConfig>();
    // sumcheck_multilinear_prod_test::<M31ExtConfig>();
}
