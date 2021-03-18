use extendr_api::prelude::*;
use peroxide::numerical::integral::*;

#[inline(always)]
fn crmh_(a: &f64, // NB: *NOT* vectorized on a
	 x: &[f64],
	 y: &[i32],
	 w: &[f64],
	 s: f64,
	 b: i32) -> f64 {
    let mut v = a.powi(b) * (-0.5*a/s).exp();
    if v.is_infinite() { return 0.0; }
    for i in 0 .. y.len() {
	let p_i = x[i].powf(a.exp()); // 'power model' CRM
	v = v * if y[i] == 0 { 1.0 - w[i] * p_i } else { p_i };
    }
    v
}

/// Integrate one of the power-model moments
///
/// @inheritParams rcrmh
/// @param b Integer in {0,1,2} telling which moment of posterior to compute
/// @export
#[extendr]
fn icrm(x: &[f64],
	y: &[i32],
	w: &[f64],
	s: f64,
	b: i32) -> f64 {
    integrate(|u| {
	let u2 = &u.powi(2);
	let a = &u/(1.0-u2); // map u in (-1,1) --> a in (-Inf,Inf)
	let da = (1.0+u2)/(1.0-u2).powi(2); // ..with change of measure
	crmh_(&a,x,y,w,s,b)*da
    }, (-1.0, 1.0), Integral::G30K61(0.000001)) // Integral::G25K51(0.0000001))
}


// Vectorize crmh1 on the 'a' parameter
fn rcrmh_(a: &[f64],
	  x: &[f64],
	  y: &[i32],
	  w: &[f64],
	  s: f64,
	  b: i32) -> Robj {
    // TODO: Assert that x, y and w all have same length?
    // TODO: Assert that y is always in {0,1}?
    let v = a.iter().map(|a| crmh_(a,x,y,w,s,b));
    v.collect_robj()
}

/// Rust implementation of \code{dfcrm::crmh*} integrands
///
/// @param a Numeric vector of evaluation points
/// @param x Numeric vector of dose-wise prior probabilities of toxicity
/// @param y Integer vector of patient-wise 0/1 toxicity indicators
/// @param w Patient-wise weights (used for TITE CRM)
/// @param s Scalar scale factor
/// @describeIn rcrmh Posterior for 1-parameter empiric (aka 'power') model
/// @export
#[extendr]
fn rcrmh(a: &[f64], x: &[f64], y: &[i32], w: &[f64], s: f64) -> Robj {
    rcrmh_(a, x, y, w, s, 0)
}

/// @describeIn rcrmh Integrand for 1st moment of empiric posterior
/// @export
#[extendr]
fn rcrmht(a: &[f64], x: &[f64], y: &[i32], w: &[f64], s: f64) -> Robj {
    rcrmh_(a, x, y, w, s, 1)
}

/// @describeIn rcrmh Integrand for 2nd moment of empiric posterior
/// @export
#[extendr]
fn rcrmht2(a: &[f64], x: &[f64], y: &[i32], w: &[f64], s: f64) -> Robj {
    rcrmh_(a, x, y, w, s, 2)
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod precautionary;
    fn icrm;
    fn rcrmh;
    fn rcrmht;
    fn rcrmht2;
}
