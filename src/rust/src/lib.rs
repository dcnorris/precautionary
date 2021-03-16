use extendr_api::prelude::*;

/// Return string `"Hello world!"` to R.
/// @export
#[extendr]
fn hello_world() -> &'static str {
    "Hello world!"
}

/// Return string `"Adios, C!"` to R.
/// @export
#[extendr]
fn adios() -> &'static str {
    "Adios, C!"
}

// A utility function not exported
fn rcrmh_(a: &[f64],
	  x: &[f64],
	  y: &[f64],
	  w: &[f64],
	  s: &[f64],
	  b: i32) -> Robj {
    // TODO: Assert that x, y and w all have same length?
    let v = a.iter().map(|a| {
	let mut v_ = a.powi(b) * (-0.5*a/s[0]).exp();
	for i in 0 .. y.len() {
	    let p_i = x[i].powf(a.exp()); // 'power model' CRM
	    v_ = v_ * if y[i] == 0.0 { 1.0 - w[i] * p_i } else { p_i };
	}
	v_
    });
    v.collect_robj()
}

/// A Rust implementation of the dfcrm::crmh posterior, which I hope will prove
/// faster to integrate() than my tuned-up (2x) reimplementation of Ken Cheung's
/// R code from package 'dfcrm'.
/// @export
#[extendr]
fn rcrmh(a: &[f64], x: &[f64], y: &[f64], w: &[f64], s: &[f64]) -> Robj {
    rcrmh_(a, x, y, w, s, 0)
}

// Posterior times x
#[extendr]
fn rcrmht(a: &[f64], x: &[f64], y: &[f64], w: &[f64], s: &[f64]) -> Robj {
    rcrmh_(a, x, y, w, s, 1)
}

// Posterior times x^2
#[extendr]
fn rcrmht2(a: &[f64], x: &[f64], y: &[f64], w: &[f64], s: &[f64]) -> Robj {
    rcrmh_(a, x, y, w, s, 2)
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod precautionary;
    fn hello_world;
    fn adios;
    fn rcrmh;
    fn rcrmht;
    fn rcrmht2;
}
