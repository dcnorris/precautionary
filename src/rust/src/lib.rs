use extendr_api::prelude::*;
use peroxide::numerical::integral::*;

#[inline(always)]
fn crmh_(a: &f64, // NB: *NOT* vectorized on a
	 x: &[f64],
	 y: &[i32],
	 w: &[f64],
	 s: f64,
	 b: i32) -> f64 {
    let mut v = a.powi(b) * (-0.5*(a/s).powi(2)).exp();
    if v.is_infinite() { return 0.0; }
    for i in 0 .. y.len() {
	let p_i = x[i].powf(a.exp()); // 'power model' CRM
	v = v * if y[i] == 0 { 1.0 - w[i] * p_i } else { p_i };
    }
    v
}

/// Rust quadrature for moments of the empiric model posterior
///
/// To match the QAGI routine used by R's \code{integrate(f, lower = -Inf, upper = Inf)},
/// the same $x = (1-t)/t$ transformation used in QUADPACK's QAGI routine is employed here,
/// albeit with 31-point Gauss-Kronrod quadrature instead of the 15-point GK reportedly
/// used in QUADPACK.
/// @seealso \url{https://en.wikipedia.org/wiki/QUADPACK#General-purpose_routines}
/// @inheritParams crmh
/// @param b Integer in {0,1,2} telling which moment of posterior to compute
/// @export
#[extendr]
fn icrm(x: &[f64],
	y: &[i32],
	w: &[f64],
	s: f64,
	b: i32) -> f64 {
    integrate(|u| {
	let a = (1.0 - &u)/&u; // map u in (0,1) --> a in (0,Inf)
	let da = &u.powi(-2); // ..with change of measure
	let fa  = crmh_(& a,x,y,w,s,b); // we integrate twice the even part of f
	let f_a = crmh_(&-a,x,y,w,s,b); // (i.e., f(a)+f(-a) over the right half
	(fa + f_a)*da
    }, (0.0, 1.0), Integral::G15K31(1e-10)) // Integral::G25K51(0.0000001))
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
    fn icrm;
    fn crmh;
    fn crmht;
    fn crmht2;
}
