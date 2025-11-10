use ark_bn254::{Bn254, Fr, G1Projective, G2Projective};
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_ff::{Field, One, Zero};
use ark_poly::{UVPolynomial, univariate::DensePolynomial};

use crate::scs::CRS;

/// We enforce that a dedicated slot w[idx_one] == 1.
/// Prover returns [Q0(τ)]_1 for (B(X) - 1) = Q0(X)*(X - D[idx_one]).
pub struct NonZeroProof {
    pub q0_tau_1: G1Projective,
    pub w_tau_2: G2Projective,  // reuse same [B(τ)]_2 commitment
}

pub fn nonzero_prove(crs: &CRS, w: &[Fr], idx_one: usize) -> NonZeroProof {
    // Build B(X) and commit
    let B = crs.interpolate(w);
    let w_tau_2 = crs.commit_poly_g2(B.coeffs());

    // KZG open at point D[idx_one] with claimed value 1: build Q0 = (B(X)-1)/(X - d)
    let d = crs.domain.element(idx_one);
    let lin = DensePolynomial::from_coefficients_vec(vec![-d, Fr::one()]);
    let mut B_minus_1 = B.clone();
    let mut c = B_minus_1.coeffs().to_vec();
    c[0] -= Fr::one();
    B_minus_1 = DensePolynomial::from_coefficients_vec(c);

    let (Q0, rem) = B_minus_1.divide_with_q_and_r(&lin).expect("div");
    debug_assert!(rem.degree() == 0 && rem.coeffs()[0].is_zero());

    let q0_tau_1 = crs.commit_poly_g1(Q0.coeffs());
    NonZeroProof { q0_tau_1, w_tau_2 }
}

// Extra GT coordinate slots for A_LV · π = b_LV:
//
// c8 = e(g1, w_tau_2)
// c9 = e(q0_tau_1, (tau - d)_2)
pub fn nonzero_verify(crs: &CRS, pi: &NonZeroProof, idx_one: usize) -> bool {
    let d = crs.domain.element(idx_one);
    // Check: [B(τ)]_2 ◦ [1]_1 = [Q0(τ)]_1 ◦ ([τ]_2 - [d]_2) + [1]_2 ◦ [1]_1
    // Rearranged linear pairing identity (proof-vs-CRS only)
    let lhs = <Bn254 as Pairing>::pairing(<Bn254 as Pairing>::G1::generator(), pi.w_tau_2);
    let term_q = <Bn254 as Pairing>::pairing(pi.q0_tau_1, crs.g2_tau_pow(1) - <Bn254 as Pairing>::G2::generator().mul_bigint(d.into_bigint()));
    let rhs = <Bn254 as Pairing>::pairing(<Bn254 as Pairing>::G1::generator(), <Bn254 as Pairing>::G2::generator()) * term_q;

    lhs == rhs
}
