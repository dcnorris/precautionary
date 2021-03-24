use extendr_api::prelude::*;

#[allow(non_snake_case)]

#[inline(always)]
fn crmh_non(exp_a_: &f64,  // NB: *NOT* vectorized
	    D: usize,      // number of doses tried so far with any non-tox
	    lnX: &[f64],   // log-skeleton at the no-tox doses
	    notox: &[usize], // numbers of no-tox obs at doses lnX
	    w: &Vec<f64>) -> f64 {
    // The non-tox factor is the nontrivial part of the CRM power model integrand,
    // because ln does not commute with any operations inside (1 - w * x^exp_a_).
    // So we must iterate at least over the weights w[.] for patients without tox.
    // BUT... this can be done without racking up per-patient transcendental ops,
    // provided we collect our factor multiplicatively (not by logs).
    let mut v_non = 1.0; // to build non-tox factor by straight multiplication
    let mut i: usize = 0; // will use this index to step thru w[i]
    for d in 0 .. D { // For each dose X[d] that was tried
	let p_d = (lnX[d]*exp_a_).exp();
	// TODO: Would an into_iter of w be the more idiomatic expression?
	//       Certainly, what I am doing is 'consuming' w, but in short
	//       bursts that are interrupted by fresh p_d calculations.
	for _ in 0 .. notox[d] {
	    v_non *= 1.0 - w[i]*p_d;
	    i += 1; // step to next weight
	}
    }

    v_non
}

// Vectorize crmh1 on the 'a' parameter
fn crmh_v(a: &[f64],
	  x: &[f64],
	  y: &[i32],
	  w: &[f64],
	  s: f64,
	  b: i32) -> Robj {
    // TODO: Assert that x, y and w all have same length?
    // TODO: Assert that y is always in {0,1}?
    let mut X: Vec<f64> = x.to_vec();
    // NB: -^ doses expressed on a prior-prob scale!
    X.sort_by(|a, b| a.partial_cmp(b).unwrap());
    X.dedup();

    const D: usize = 10; // for statically allocating several arrays

    let mut lnX_: [f64; 10] = [0.0; D];
    for d in 0 .. X.len() {
	lnX_[d] = X[d].ln();
    }

    // It also proves useful to tally dose-wise sums of toxicities:
    let mut Y: [f64; 10] = [0.0; D];
    for i in 0 .. y.len() {
	match X.iter().position(|&p| p==x[i]) {
	    Some(d) => {
		if y[i]==1 {
		    Y[d] += 1.0; // tally toxicities
		}
	    }
	    None => {}
	}
    }
    // NB: Although I could eliminate Y, by computing v_tox inside the above loop,
    //     my aim in introducing Y is partly to prepare the way for eventually
    //     accepting it as a parameter to this function.

    // We can calculate the toxic term irrespective of 'a' (except for its exp(a) factor)
    let mut log_vtox = 0.0;
    for d in 0 .. D { // TODO: Use an iterator-based expression?
	log_vtox += Y[d]*lnX_[d];
    }

    // To speed up the no-tox factor iteration in crmh_, we now precompute
    // some details of that loop.
    // Initially, let me treat the w==1 case, where all I need is a D-vector of
    // counts of no-tox cases.
    let mut notox: [usize; 10] = [0; D];
    let mut W: Vec<f64> = w.to_vec(); // TODO: Rationalize the sizing
    let mut iv: i32 = 0; // just a temp var to make a check
    let mut iw: usize = 0;
    for d in 0 .. X.len() { // We make a single pass through unique doses tried
	for i in 0 .. y.len() { // ..looking for patients
	    if y[i] == 0 && x[i] == X[d] { // ..who got that dose, without toxicity.
		notox[d] += 1; // For each such patient we increment a dose-wise counter,
		W[iw] = w[i];  // store that patient's weight in vector W.
		iw += 1;       // and point to the next empty spot in W.
		iv += 1; // (temporarily track iw as an i32 for comparison)
	    }
	}
    }

    let v = a.iter().map(|a| if a > &709.0 { 0.0 } else {
	let log_vconst = -0.5 * (a/s).powi(2);
	let exp_a_ = a.exp();
	crmh_non(&exp_a_,X.len(),&lnX_,&notox,&W)
	    * (log_vconst + exp_a_*log_vtox).exp()
	    * a.powi(b)
    });
    v.collect_robj()
}

/// Rust implementation of \code{dfcrm::crmh*} integrands
///
/// @param a Numeric vector of evaluation points
/// @param x Numeric vector of dose-wise prior probabilities of toxicity
/// @param y Integer vector of patient-wise 0/1 toxicity indicators
/// @param w Patient-wise weights (used for TITE CRM)
/// @param s Scalar scale factor
/// @describeIn crmh Posterior for 1-parameter empiric (aka 'power') model
/// @export
#[extendr]
fn crmh(a: &[f64], x: &[f64], y: &[i32], w: &[f64], s: f64) -> Robj {
    crmh_v(a, x, y, w, s, 0)
}

/// @describeIn crmh Integrand for 1st moment of empiric posterior
/// @export
#[extendr]
fn crmht(a: &[f64], x: &[f64], y: &[i32], w: &[f64], s: f64) -> Robj {
    crmh_v(a, x, y, w, s, 1)
}

/// @describeIn crmh Integrand for 2nd moment of empiric posterior
/// @export
#[extendr]
fn crmht2(a: &[f64], x: &[f64], y: &[i32], w: &[f64], s: f64) -> Robj {
    crmh_v(a, x, y, w, s, 2)
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod precautionary;
    fn crmh;
    fn crmht;
    fn crmht2;
}
