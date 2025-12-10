//src/iip.rs
use ark_bn254::{Bn254, Fr, G1Projective, G2Projective};
use ark_ec::PrimeGroup;
use ark_ec::pairing::Pairing;
use ark_ff::{Field, One, PrimeField, Zero};
use ark_poly::{DenseUVPolynomial, Polynomial, univariate::DensePolynomial};
use crate::helpers::{add_constant, sub_poly, scale_poly, mul_by_xk, mul_poly, poly_from_coeffs, div_rem};

use crate::scs::CRS;

/// Public digest (vk) for IIP, as in Construction 6.
#[allow(non_snake_case)]
#[allow(dead_code)]
#[derive(Clone)]
pub struct IIPDigest {
    pub x_star: Fr,                           // we use 0
    pub y_star: Fr,                           // 1/n
    pub C: G1Projective,                      // y* · [Σ s_i L_i(τ)]_1
    pub Z_tau_2: G2Projective,                // [Z(τ)]_2
    pub tau_2: G2Projective,                  // [τ]_2  (since x*=0, [τ - x*]_2 = [τ]_2)
    pub tau_N_minus_n_plus_2_2: G2Projective, // [τ^{N-n+2}]_2
    pub tau_N_2: G2Projective,                // [τ^N]_2
    pub n: usize,
    pub N: usize,
}

#[derive(Clone)]
#[allow(non_snake_case)]
pub struct IIPProof {
    pub w_tau_2: G2Projective,      // [B(τ)]_2 = SCS(G2).Commit(w)
    pub v_g1: G1Projective,         // v = Σ w_i [s_i]_1
    pub QZ_tau_1: G1Projective,     // [Q_Z(τ)]_1
    pub QX_tau_1: G1Projective,     // [Q_X(τ)]_1
    pub QX_hat_tau_1: G1Projective, // [Q̂_X(τ)]_1 = [X^{N-n+2} Q_X(X)]_1
    pub v_hat_tau_1: G1Projective,  // [v̂(τ)]_1 = [X^N · (Σ w_i s_i)]_1
}

/// Build vk for IIP given public index s in F^n
#[allow(non_snake_case)]
pub fn iip_digest(crs: &CRS, s: &[Fr]) -> IIPDigest {
    assert_eq!(s.len(), crs.n);
    // A(X) interpolates s over D
    let A = crs.interpolate(s);
    let A_coeffs = A.coeffs();

    
    let A_tau_1 = crs.commit_poly_g1(A_coeffs);
    //let C = A_tau_1.mul_bigint(crs.n_inv.into_bigint()); // this is the paper’s scaled variant we (Construction 6), we must refactor so that the pairing identity balances well.
    let C = A_tau_1;
    

    // [Z(τ)]_2 from precomputed coeffs
    let Z_tau_2 = crs.commit_poly_g2(&crs.vanishing_coeffs);

    IIPDigest {
        x_star: Fr::zero(),
        y_star: crs.n_inv,
        C,
        Z_tau_2,
        tau_2: crs.g2_tau_pow(1),
        tau_N_minus_n_plus_2_2: crs.g2_tau_pow(crs.N - crs.n + 1),
        //tau_N_minus_n_plus_2_2: crs.g2_tau_pow(crs.N - crs.n + 2),
        tau_N_2: crs.g2_tau_pow(crs.N),
        n: crs.n,
        N: crs.N,
    }
}

/// Prover: compute B(X), v, Q_X, Q_Z, and the “hatted” terms.
#[allow(non_snake_case)]
pub fn iip_prove(crs: &CRS, s: &[Fr], w: &[Fr]) -> IIPProof {
    assert_eq!(s.len(), crs.n);
    assert_eq!(w.len(), crs.n);

    // A(X), B(X)
    let A = crs.interpolate(s);
    let B = crs.interpolate(w);

    // Commit w
    let w_tau_2 = crs.commit_poly_g2(B.coeffs());

    // v = Σ w_i [s_i]_1
    let mut v_scalar = Fr::zero();
    for (wi, si) in w.iter().zip(s.iter()) {
        v_scalar += *wi * *si;
    }
    let v_g1 = <Bn254 as Pairing>::G1::generator().mul_bigint(v_scalar.into_bigint());

    // P(X) = A(X)B(X) - (Σ w_i s_i)/y*
    let mut P = mul_poly(&A, &B);
    //let t = v_scalar * crs.n_inv.inverse().unwrap();  
    //let t = v_scalar * crs.n_inv;  
    let n_field = crs.n_inv.inverse().unwrap(); 
    let t = v_scalar * n_field; 
    // subtract constant t
    let mut P_coeffs = P.coeffs().to_vec();
    if P_coeffs.is_empty() {
        P_coeffs.push(-t);
    } else {
        P_coeffs[0] -= t;
    }
    P = poly_from_coeffs(P_coeffs);

    // Z(X)
    let Z = DensePolynomial::from_coefficients_vec(crs.vanishing_coeffs.clone());

    // 1) Divide P by Z: P = QZ * Z + R, deg R < n
    let (mut QZ, mut R) = div_rem(&P, &Z);

    // 2) Adjust so that R(x*) = 0 with x* = 0:
    let x_star = Fr::zero();
    let Z_x = Z.evaluate(&x_star);
    let R_x = R.evaluate(&x_star);
    if !R_x.is_zero() {
        // set QZ' = QZ + c with c = R(0)/Z(0), R' = R - c*Z ⇒ R'(0)=0
        let c = R_x * Z_x.inverse().unwrap();
        QZ = add_constant(&QZ, c);
        R = sub_poly(&R, &scale_poly(&Z, c));
    }

    // 3) Now R is divisible by (X - x*), define QX = R / (X - x*)
    let lin = DensePolynomial::from_coefficients_vec(vec![-x_star, Fr::one()]);
    let (QX, rem) = div_rem(&R, &lin);
    debug_assert!(rem.is_zero(), "R(X) not divisible by (X - x*)");

    // Hatted polynomials:
    // Q̂_X(X) = X^{N-n+1} Q_X(X)
    let QX_hat = mul_by_xk(&QX, (crs.N - crs.n + 1) as usize);
    //let QX_hat = mul_by_xk(&QX, (crs.N - crs.n + 2) as usize);
    
    // v̂(X) = X^N * (Σ w_i s_i)  (a pure monomial with that coefficient)
    let mut vhat_coeffs = vec![Fr::zero(); crs.N + 1];
    vhat_coeffs[crs.N] = v_scalar;
    let vhat = DensePolynomial::from_coefficients_vec(vhat_coeffs);

    IIPProof {
        w_tau_2,
        v_g1,
        QZ_tau_1: crs.commit_poly_g1(QZ.coeffs()),
        QX_tau_1: crs.commit_poly_g1(QX.coeffs()),
        QX_hat_tau_1: crs.commit_poly_g1(QX_hat.coeffs()),
        v_hat_tau_1: crs.commit_poly_g1(vhat.coeffs()),
    }
}

/// Verifier: the three linear checks from Construction 6 (no proof–proof pairing).
/// // GT coordinate slots used for A_LV · π = b_LV:
//
// c0 = e(C, w_tau_2)
// c1 = e(v_g1, y_star^{-1} * g2)
// c2 = e(QX_tau_1, tau_2)
// c3 = e(QZ_tau_1, Z_tau_2)
// c4 = e(QX_tau_1, tau_N_minus_n_plus_2_2)
// c5 = e(QX_hat_tau_1, g2)
// c6 = e(v_g1, tau_N_2)
// c7 = e(v_hat_tau_1, g2)
// (NonZero adds c8,c9 in nonzero.rs)
#[allow(non_snake_case)]
pub fn iip_verify(d: &IIPDigest, pi: &IIPProof) -> bool {
    // 1) C ◦ w = v ◦ [y*^{-1}]_2 + [QX(τ)]_1 ◦ [τ - x*]_2 + [QZ(τ)]_1 ◦ Z
    let lhs1 = <Bn254 as Pairing>::pairing(d.C, pi.w_tau_2);

    // v ◦ [y*^{-1}]_2
    let y_inv = d.y_star.inverse().unwrap();
    let v_g1_scaled = pi.v_g1.mul_bigint(y_inv.into_bigint());
    let rhs1_v = <Bn254 as Pairing>::pairing(v_g1_scaled, <Bn254 as Pairing>::G2::generator());

    // [QX(τ)]_1 ◦ [τ - x*]_2, and we have x* = 0 ⇒ [τ - x*]_2 = [τ]_2
    let term_qx = <Bn254 as Pairing>::pairing(pi.QX_tau_1, d.tau_2);

    // [QZ(τ)]_1 ◦ Z
    let term_qz = <Bn254 as Pairing>::pairing(pi.QZ_tau_1, d.Z_tau_2);

    // Multiply underlying GT elements (Fq12) and wrap back into PairingOutput
    let rhs1_total = rhs1_v + term_qx + term_qz;

    if lhs1 != rhs1_total {
        return false;
    }

    // 2) [QX(τ)]_1 ◦ [τ^{N-n+2}]_2 = [Q̂X(τ)]_1 ◦ [1]_2
    let lhs2 = <Bn254 as Pairing>::pairing(pi.QX_tau_1, d.tau_N_minus_n_plus_2_2);
    let rhs2 = <Bn254 as Pairing>::pairing(pi.QX_hat_tau_1, <Bn254 as Pairing>::G2::generator());
    if lhs2 != rhs2 {
        return false;
    }

    // 3) v ◦ [τ^N]_2 = [v̂(τ)]_1 ◦ [1]_2
    let lhs3 = <Bn254 as Pairing>::pairing(pi.v_g1, d.tau_N_2);
    let rhs3 = <Bn254 as Pairing>::pairing(pi.v_hat_tau_1, <Bn254 as Pairing>::G2::generator());
    if lhs3 != rhs3 {
        return false;
    }

    true
}