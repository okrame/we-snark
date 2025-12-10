// src/helpers.rs

use ark_bn254::Fr;
use ark_ff::{Zero};
use std::ops::Mul;
use ark_poly::{DenseUVPolynomial, univariate::DensePolynomial};

/// Add a constant to a polynomial: p(X) + c
pub fn add_constant(p: &DensePolynomial<Fr>, c: Fr) -> DensePolynomial<Fr> {
    let mut v = p.coeffs().to_vec();
    if v.is_empty() {
        v.push(c);
    } else {
        v[0] += c;
    }
    DensePolynomial::from_coefficients_vec(v)
}

/// Subtract two polynomials: a(X) - b(X)
pub fn sub_poly(a: &DensePolynomial<Fr>, b: &DensePolynomial<Fr>) -> DensePolynomial<Fr> {
    DensePolynomial::from_coefficients_vec((a - b).coeffs().to_vec())
}

/// Scale a polynomial by a constant: c * p(X)
pub fn scale_poly(p: &DensePolynomial<Fr>, c: Fr) -> DensePolynomial<Fr> {
    DensePolynomial::from_coefficients_vec(
        p.coeffs().iter().map(|x| *x * c).collect()
    )
}

/// Multiply polynomial by X^k: X^k * p(X)
pub fn mul_by_xk(p: &DensePolynomial<Fr>, k: usize) -> DensePolynomial<Fr> {
    let mut v = vec![Fr::zero(); k];
    v.extend_from_slice(p.coeffs());
    DensePolynomial::from_coefficients_vec(v)
}

/// Multiply two polynomials
pub fn mul_poly(a: &DensePolynomial<Fr>, b: &DensePolynomial<Fr>) -> DensePolynomial<Fr> {
        a.mul(b)
    }

/// Create polynomial from coefficient vector
pub fn poly_from_coeffs(coeffs: Vec<Fr>) -> DensePolynomial<Fr> {
        DensePolynomial::from_coefficients_vec(coeffs)
}

/// Polynomial division with remainder: returns (quotient, remainder)
/// where dividend = quotient * divisor + remainder
#[allow(non_snake_case)]
pub fn div_rem(
        P: &DensePolynomial<Fr>,
        Q: &DensePolynomial<Fr>,
    ) -> (DensePolynomial<Fr>, DensePolynomial<Fr>) {
        let q = P / Q;
        let r = P - &(&q * Q);
        (q, r)
    }