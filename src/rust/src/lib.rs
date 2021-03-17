use extendr_api::prelude::*;
use peroxide::numerical::integral::*;

/// Return string `"Adios, C!"` to R.
/// @export
#[extendr]
fn adios() -> &'static str {
    "Adios, C!"
}

#[inline(always)]
fn crmh_(a: &f64, // NB: *NOT* vectorized on a
	 x: &[f64],
	 y: &[i32],
	 w: &[f64],
	 s: f64,
	 b: i32) -> f64 {
    let mut v = a.powi(b) * (-0.5*a/s).exp();
    for i in 0 .. y.len() {
	let p_i = x[i].powf(a.exp()); // 'power model' CRM
	v = v * if y[i] == 0 { 1.0 - w[i] * p_i } else { p_i };
    }
    v
}

/// Integrate one of the power-model moments
/// @export
#[extendr]
fn icrm(x: &[f64],
	y: &[i32],
	w: &[f64],
	s: f64,
	b: i32) -> f64 {
    integrate(|u| crmh_(&u,x,y,w,s,b), (-10.0, 10.0),
	      Integral::G7K15(0.000001))
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

/// A Rust implementation of the dfcrm::crmh posterior, which I hope will prove
/// faster to integrate() than my tuned-up (2x) reimplementation of Ken Cheung's
/// R code from package 'dfcrm'.
/// @export
#[extendr]
fn rcrmh(a: &[f64], x: &[f64], y: &[i32], w: &[f64], s: f64) -> Robj {
    rcrmh_(a, x, y, w, s, 0)
}

/// Posterior times x
/// @export
#[extendr]
fn rcrmht(a: &[f64], x: &[f64], y: &[i32], w: &[f64], s: f64) -> Robj {
    rcrmh_(a, x, y, w, s, 1)
}

/// Posterior times x^2
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
    fn adios;
    fn icrm;
    fn rcrmh;
    fn rcrmht;
    fn rcrmht2;
}
