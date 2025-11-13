//src/nonzero.rs
use ark_bn254::{Bn254, Fr, G1Projective, G2Projective};
use ark_ec::pairing::Pairing;
use ark_ec::PrimeGroup;
use ark_ff::{One, PrimeField, Zero};
use ark_poly::{
    univariate::DensePolynomial,
    DenseUVPolynomial,
    EvaluationDomain
};

use crate::scs::CRS;

/// We enforce that a dedicated slot w[idx_one] == 1.
/// Prover returns [Q0(τ)]_1 for (B(X) - 1) = Q0(X)*(X - D[idx_one]).
#[derive(Clone)]
pub struct NonZeroProof {
    pub q0_tau_1: G1Projective,
    pub w_tau_2: G2Projective, // reuse same [B(τ)]_2 commitment
}

/// Synthetic division by (X - d) for polynomials with coefficients from lowest to highest degree.
/// Given P(X) = sum_i c[i] X^i, returns (Q(X), r) such that:
/// P(X) = (X - d) Q(X) + r
fn divide_by_linear(poly: &DensePolynomial<Fr>, d: Fr) -> (DensePolynomial<Fr>, Fr) {
    let coeffs = poly.coeffs();
    let n = coeffs.len();

    if n == 0 {
        return (DensePolynomial::zero(), Fr::zero());
    }
    if n == 1 {
        // constant polynomial
        return (DensePolynomial::zero(), coeffs[0]);
    }

    // Convert to descending coefficients a[0]..a[n] (degree n .. 0)
    let a: Vec<Fr> = coeffs.iter().cloned().rev().collect();

    // Standard synthetic division in descending convention
    let mut b = vec![Fr::zero(); n];
    b[0] = a[0];
    for i in 1..n {
        // b_i = a_i + d * b_{i-1}
        b[i] = a[i] + d * b[i - 1];
    }
    let r = b[n - 1]; // remainder
    let quot_desc = &b[..n - 1];

    // Convert quotient back to ascending order
    let q_coeffs: Vec<Fr> = quot_desc.iter().cloned().rev().collect();
    let q = DensePolynomial::from_coefficients_vec(q_coeffs);

    (q, r)
}

#[allow(non_snake_case)]
pub fn nonzero_prove(crs: &CRS, w: &[Fr], idx_one: usize) -> NonZeroProof {
    // Build B(X) and commit
    let B = crs.interpolate(w);
    let w_tau_2 = crs.commit_poly_g2(B.coeffs());

    // KZG open at point D[idx_one] with claimed value 1:
    // build Q0 = (B(X) - 1)/(X - d)
    let d = crs.domain.element(idx_one);

    // B_minus_1(X) = B(X) - 1
    let mut c = B.coeffs().to_vec();
    if c.is_empty() {
        c.push(-Fr::one());
    } else {
        c[0] -= Fr::one();
    }
    let B_minus_1 = DensePolynomial::from_coefficients_vec(c);

    let (Q0, rem) = divide_by_linear(&B_minus_1, d);
    debug_assert!(rem.is_zero(), "B(X) - 1 not divisible by (X - d)");

    let q0_tau_1 = crs.commit_poly_g1(Q0.coeffs());
    NonZeroProof { q0_tau_1, w_tau_2 }
}

// Extra GT coordinate slots for A_LV · π = b_LV:
//
// c8 = e(g1, w_tau_2)
// c9 = e(q0_tau_1, (tau - d)_2)
pub fn nonzero_verify(crs: &CRS, pi: &NonZeroProof, idx_one: usize) -> bool {
    let d = crs.domain.element(idx_one);

    // [τ]_2 - [d]_2
    let tau_2 = crs.g2_tau_pow(1);
    let d_g2 =
        <Bn254 as Pairing>::G2::generator().mul_bigint(d.into_bigint());
    let tau_minus_d_2 = tau_2 - d_g2;

    // Check (in additive GT notation):
    // e(g1, [B(τ)]_2) = e(g1, [1]_2) + e([Q0(τ)]_1, [τ - d]_2)
    //
    // i.e. B(d) = 1 enforced via KZG opening
    let lhs =
        <Bn254 as Pairing>::pairing(<Bn254 as Pairing>::G1::generator(), pi.w_tau_2);
    let term_q =
        <Bn254 as Pairing>::pairing(pi.q0_tau_1, tau_minus_d_2);
    let base =
        <Bn254 as Pairing>::pairing(
            <Bn254 as Pairing>::G1::generator(),
            <Bn254 as Pairing>::G2::generator(),
        );

    // GT is modelled additively: product of pairings becomes sum in PairingOutput.
    let rhs = base + term_q;

    lhs == rhs
}
