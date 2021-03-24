use extendr_api::prelude::*;

#[allow(non_snake_case)]

#[inline(always)]
fn crmh_(a: &f64, // NB: *NOT* vectorized on a
	 lnX: &[f64],
	 D: usize,
	 Y: &[f64],
	 notox: &[i32],
	 w: &[f64],
	 s: f64,
	 b: i32) -> f64 {
    // Below was the original, very short computation.
    // I count N+1 transcendental ops, compared with 2D+1 in new implementation.
    // If an average of 4 patients are enrolled per dose level, this translates
    // to a savings of (4D - 2D)/(4D+1) ~ 50% of effort.
    // TODO: By precalculating the X[d].ln()'s in the calling function, the cost
    //       of these txops amortizes to zero! Thus we can bring this routine
    //       down to just D+1 transcendendal ops, saving 75% of effort!
    //       I might even dare to call these ln()'s the EXOSKELETON.
    /*
    let mut v = a.powi(b) * (-0.5*(a/s).powi(2)).exp();         // 1 exp()
    if v.is_infinite() { return 0.0; }
    for i in 0 .. y.len() {
	let p_i = x[i].powf(a.exp()); // 'power model' CRM      // N powf()'s
	v = v * if y[i] == 0 { 1.0 - w[i] * p_i } else { p_i };
    }
    v
     */

    // I'm going to undertake some refactoring toward an algorithm
    // that reduces transcendental ops to a bare minimum. The key
    // opportunity in this regard lies in COMPUTING x[i]^exp_a_
    // JUST ONCE PER DOSE LEVEL, rather than once per patient.
    // Some bookkeeping preliminaries are needed to make this possible.
    // As part of my initial refactoring, I will do this initially here
    // in this routine; but this computation is one that ought to be
    // bubbled upward as far as possible.

    if a > &709.0 { return 0.0 } // "Not gonna exp() it; wouldn't be prudent."
    let exp_a_ = a.exp();                                             // 1 exp()

    // The objective function can factored as: v = vconst * log_vtox.exp() * v_non.
    // (Let's call log_vtox the 'toxic term' and v_non the 'non-tox factor'.)
    let log_vconst = -0.5 * (a/s).powi(2);
    if log_vconst > 709.0 { return 0.0; } // saves time, avoids returning Inf*0=NaN

    let mut log_vtox = 0.0;
    for d in 0 .. D { // TODO: Use an iterator-based expression?
	log_vtox += Y[d]*lnX[d];
    }
    log_vtox = log_vtox * exp_a_; // this completes the toxic term

    // The non-tox factor is more difficult, because ln does not
    // commute with any operations inside (1 - w * x^exp_a_) 8^(.
    // So we must iterate over the patients. BUT... this can be
    // done without racking up per-patient transcendental ops,
    // provided we collect our factor multiplicatively (not by logs).
    // The key point is that the following computation costs only
    // X.len() transcendental ops -- just 1 powf() per dose tried.
    let mut v_non = 1.0; // to build non-tox factor by straight multiplication
    for d in 0 .. D { // For each dose X[d] that was tried
	let p_d = (lnX[d]*exp_a_).exp();                              // D powf()'s
	for i in 0 .. notox[d] {
	    v_non *= (1.0 - 1.0*p_d) // TODO: support general case as below:
	    //v_non *= (1.0 - w[i]*p_d)
	}
    }

    let vfast = a.powi(b) * (log_vconst + log_vtox).exp() * v_non;    // 1 exp()
    vfast
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

    let mut X_: [f64; 10] = [0.0; D];
    let mut lnX_: [f64; 10] = [0.0; D];
    for d in 0 .. X.len() {
	X_[d] = X[d];
	lnX_[d] = X_[d].ln();
    }

    // It also proves useful to tally dose-wise sums of toxicities:
    let mut Y: [f64; 10] = [0.0; 10]; // Static allocation sufficient for 10 doses
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

    // To speed up the no-tox factor iteration in crmh_ (also eliminating the need
    // for passing X_ at all!), we now precompute some details of that loop.
    // Initially, let me treat the w==1 case, where all I need is a D-vector of
    // counts of no-tox cases.
    let mut notox: [i32; 10] = [0; D];
    for d in 0 .. X.len() {
	for i in 0 .. y.len() {
	    if y[i] == 0 && x[i] == X[d] {
		notox[d] += 1;
	    }
	}
    }

    let v = a.iter().map(|a| crmh_(a,&lnX_,X.len(),&Y,&notox,w,s,b));
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
