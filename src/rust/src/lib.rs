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

// Decode an f64-encoded trial-path tally to dose-wise numbers of tox & no-tox.
// TODO: Generalize to cohorts with n=1 or n=2, ideally by encoding n itself
//       within the 'obs'.
#[inline(always)]
fn decode_cohorts(obs: f64,
		  tox: &mut [i32; 8], // dose-wise toxicity tally
		  nos: &mut [i32; 8]  // dose-wise no-tox tally
) -> () {
    let iobs = obs as i32;  // integer part of obs
    let fobs = obs.fract(); // fractional part
    let mut d6 = 1;
    let mut d_16 = 1.0;
    for d in 0 .. 8 {
	let enr_d = (iobs.rem_euclid(6*d6)).div_euclid(d6); // d6 is 6^d
	d6 = 6*d6;
	// Extract tox counts from fractional part of 'obs'
	d_16 = d_16 / 16.0; // d_16 is 16^-(d+1)
	tox[d] = fobs.rem_euclid(16.0*d_16).div_euclid(d_16) as i32;
	// Compute non-toxicities
	nos[d] = 3*enr_d - tox[d]; // NB: n=3 is hard-coded here
    }
}

// TODO: Try mul_add for cleaner implementation.
// TODO: Consider bitwise construction using f64::from_bits.
fn encode_cohorts(enr: &[i32], // dose-wise number of cohorts enrolled
		  tox: &[i32]  // dose-wise number of toxicities
) -> f64 {
    let mut obs = 0.0;
    let mut eix = 1.0;
    let mut tix = 1.0/16.0;
    for d in 0 .. 8 {
	obs += eix * enr[d] as f64;
	eix = eix * 6.0;
	obs += tix * tox[d] as f64;
	tix = tix / 16.0;
    }
    obs
}

// A [w==1] version of crmh_v that takes f64-encoded trial observations as parameter.
// When all w's are 1.0, the patients become exchangeable, and only *counts* (which
// are encoded in the f64 'obs') need to be supplied.
// (This could also be seen from the converse perspective: the non-exchangeability
// of patients 'breaks' the encoding used here, necessitating the more complicated
// routine implmented in crmh_v.)
// The encoding assumes a skeleton of 8 doses at most, which I allocate statically.
// Note that this implementation allows the caller to provide the raw log-skeleton
// without any path-specific pre-processing, which is a benefit!
// One way to see this function in general, is as exploring a special case of the
// general concept of creating atomic indexes for caching expensive computations.
/// Rust implementation of \code{dfcrm::crmh*} integrands for w==1 case
///
/// @param a Numeric vector of evaluation points
/// @param obs An f64-encoding of dose-wise enrollment and toxicity counts
/// @param ln_x A numeric vector of dose-wise prior log-probabilities of toxicity
/// @param s Scalar scale factor
/// @param b Order of moment to calculate (0, 1 or 2)
/// @export
#[extendr]
fn crmh_ev(a: &[f64],
	   obs: f64, // encoding of dose-wise enrolled cohorts & toxicity counts
	   ln_x: &[f64], // DOSE-WISE log-skeleton -- different from crmh_v!
	   s: f64,
	   b: i32) -> Robj {

    const D: usize = 8;
    let mut enr: [i32; D] = [0; D]; // enrollment (in cohorts) at each dose
    let tox: &mut [i32; D] = &mut [0; D]; // toxicities at each dose
    let nos: &mut [i32; D] = &mut [0; D]; // non-tox count at each dose

    decode_cohorts(obs, tox, nos);

    // As a check, reconstitute the 'obs'
    for d in 0 .. 8 {
	enr[d] = (tox[d] + nos[d])/3;
    }
    let check = encode_cohorts(&enr, tox);
    let diff = (obs - check).abs();
    assert!(diff < 1e-15, "Ouch! obs = {}, check = {}", obs, check);

    // As in crmh_v, this calculation is something that could be done
    // by the caller. But note that, with ln_x being here DOSE-INDEXED,
    // the proper approach to the calculation is somewhat different.
    // So I demonstrate that approach here for the caller's benefit.
    let mut log_vtox = 0.0;
    for d in 0 .. ln_x.len() {
	log_vtox +=  ln_x[d] * tox[d] as f64; // TODO: mul_add
    }

    let v = a.iter().map(|a| if a > &709.0 { 0.0 } else {
	let log_vconst = -0.5 * (a/s).powi(2);
	let exp_a_ = a.exp();
	let o_factor = { // x/o-notation-inspired name for 'non-tox factor'
	    let mut nontox_ = 1.0;
	    for d in 0 .. ln_x.len() {
		if nos[d] > 0 { // skip over any 100%-toxic dose
		    let p_d = (ln_x[d]*exp_a_).exp(); // TXOP
		    nontox_ *= (1.0 - p_d).powi(nos[d]);
		}
	    }
	    nontox_
	};
	let x_term = log_vconst + exp_a_*log_vtox;
	assert!(o_factor.is_finite(), "O-factor = {}; obs = {}, a = {}", o_factor, obs, a);
	assert!(x_term.exp().is_finite(), "X-term = {}", x_term);
	o_factor * x_term.exp() * a.powi(b)
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
    fn crmh_ev;
    fn crmh;
    fn crmht;
    fn crmht2;
}
