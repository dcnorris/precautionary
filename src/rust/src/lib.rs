use extendr_api::prelude::*;

#[allow(non_snake_case)]

#[inline(always)]
fn crmh_non(exp_a_: &f64,  // NB: *NOT* vectorized
	    ln_x: &[f64],   // the log-skeleton, patient-wise
	    w: &[f64]) -> f64 {
    // The non-tox factor is the nontrivial part of the CRM power model integrand,
    // because ln does not commute with any operations inside (1 - w * x^exp_a_).
    // By encoding y in w (via w==0.0 iff y==1), however, we can make the scan
    // straightforward.
    // It is also possible to skip some exp()'s in this loop, by sorting patients
    // (i.e., ln_x's and w's) before entry, to make the ln_x's contiguous.
    let mut ln_x_ = ln_x[0];
    let mut p_ = (ln_x_*exp_a_).exp();
    let mut v_non = 1.0 - w[0]*p_;
    for i in 1 .. w.len() {
	if w[i] > 0.0 { // This 'guard' is not necessary for correctness ...
	    if ln_x[i] != ln_x_ {
		ln_x_ = ln_x[i];
		p_ = (ln_x_*exp_a_).exp(); // TXOP
	    } // (above recalc gets skipped often if ln_x contiguous)
	    v_non *= 1.0 - w[i]*p_; // ... because w[i]==0 case gives *= 1.0 here.
	} // This guard's performance impacts are TBD, and likely data-dependent.
    }

    v_non
}

// Vectorize crmh1 on the 'a' parameter
fn crmh_v(a: &[f64],
	  ln_x: &[f64], // patient-wise log-skeleton
	  w: &[f64], // patient-wise weights; tox iff 0.0
	  s: f64,
	  b: i32) -> Robj {

    // We can calculate the toxic term irrespective of 'a',
    // except for its exp(a) factor which we multiply on at the end.
    let mut log_vtox = 0.0;
    for i in 0 .. w.len() {
	if w[i] == 0.0 { // iterate over toxicities
	    log_vtox += ln_x[i];
	}
    }

    let v = a.iter().map(|a| if a > &709.0 { 0.0 } else {
	let log_vconst = -0.5 * (a/s).powi(2);
	let exp_a_ = a.exp();
	crmh_non(&exp_a_,&ln_x,&w)
	    * (log_vconst + exp_a_*log_vtox).exp()
	    * a.powi(b)
    });

    v.collect_robj()
}

/// Rust implementation of \code{dfcrm::crmh*} integrands
///
/// @param a Numeric vector of evaluation points
/// @param ln_x A numeric vector of dose-wise prior log-probabilities of toxicity
/// @param w Patient-wise weights (used for TITE CRM), also encoding toxicity
/// by \code{w[i] == 0.0}.
/// @param s Scalar scale factor
/// @describeIn crmh Posterior for 1-parameter empiric (aka 'power') model
/// @export
#[extendr]
fn crmh(a: &[f64], ln_x: &[f64], w: &[f64], s: f64) -> Robj {
    crmh_v(a, ln_x, w, s, 0)
}

/// @describeIn crmh Integrand for 1st moment of empiric posterior
/// @export
#[extendr]
fn crmht(a: &[f64], ln_x: &[f64], w: &[f64], s: f64) -> Robj {
    crmh_v(a, ln_x, w, s, 1)
}

/// @describeIn crmh Integrand for 2nd moment of empiric posterior
/// @export
#[extendr]
fn crmht2(a: &[f64], ln_x: &[f64], w: &[f64], s: f64) -> Robj {
    crmh_v(a, ln_x, w, s, 2)
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
