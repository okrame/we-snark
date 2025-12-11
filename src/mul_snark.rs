// src/mul_snark.rs

use ark_bn254::{Bn254, Fr};
use ark_ec::PrimeGroup;
use ark_ec::pairing::Pairing;
use ark_ff::{One, Zero};
use ark_poly::{DenseUVPolynomial, Polynomial, univariate::DensePolynomial};

use crate::compiler::{CompiledQAP, LVSystemParams};
use crate::gadgets::arithmetic::QAPCommit;
use crate::helpers::mul_by_xk;
use crate::iip::{iip_digest, iip_prove};
use crate::nonzero::nonzero_prove;
use crate::scs::CRS;
use crate::verifier::{LVDigest, LVProof};

/// Fixed-size MulCircuit witness: w = [x, y, z, 1].
#[derive(Clone, Debug)]
pub struct MulWitness {
    pub x: Fr,
    pub y: Fr,
    pub z: Fr,
}

impl MulWitness {
    /// Convert to the evaluation vector [x, y, z, 1] on D.
    pub fn to_vec(&self) -> Vec<Fr> {
        vec![self.x, self.y, self.z, Fr::from(1u32)]
    }
}

/// Public parameters (vk) for the LV-SNARK.
/// For now this is just a wrapper around LVDigest + the public index s.
#[derive(Clone)]
pub struct MulDigest {
    pub lv: LVDigest,
    // Selectors for the three witness slots
    pub s_x: Vec<Fr>, // [1,0,0,0]
    pub s_y: Vec<Fr>, // [0,1,0,0]
    pub s_z: Vec<Fr>, // [0,0,1,0]
}

/// Proof object for MulCircuit: reuses LVProof as-is.
#[derive(Clone)]
pub struct MulProof {
    pub lv: LVProof,
}

#[allow(non_snake_case)]
impl MulDigest {
    pub fn setup(crs: &CRS) -> Self {
        assert_eq!(
            crs.n, 4,
            "MulCircuit is currently hard-coded for n=4 (slots [x,y,z,1])"
        );

        // Selectors for x, y, z in w = [x, y, z, 1]
        let s_x = vec![
            Fr::from(1u32),
            Fr::from(0u32),
            Fr::from(0u32),
            Fr::from(0u32),
        ];
        let s_y = vec![
            Fr::from(0u32),
            Fr::from(1u32),
            Fr::from(0u32),
            Fr::from(0u32),
        ];
        let s_z = vec![
            Fr::from(0u32),
            Fr::from(0u32),
            Fr::from(1u32),
            Fr::from(0u32),
        ];

        // Z(X) = X - 1 (Mul QAP vanishing poly on the single gate)
        let z_poly = DensePolynomial::from_coefficients_vec(vec![-Fr::one(), Fr::one()]);
        let mul_z_tau_2 = crs.commit_poly_g2(z_poly.coeffs());

        // IIP vk’s for x, y, z
        let iip_vk_x = iip_digest(crs, &s_x);
        let iip_vk_y = iip_digest(crs, &s_y);
        let iip_vk_z = iip_digest(crs, &s_z);

        // Max degree bound for the SCS witness polynomial B(X) for w=[x,y,z,1]
        //let tau_N_minus_d_1 = crs._g1_tau_pow(N - d_bound);
        let sys = LVSystemParams::for_mul_demo(crs);
        let d_bound = sys.d_bound; // with n=4, d_bound=3
        let N = crs.N;
        // [τ^{N-d}]_1 in G1
        let tau_N_minus_d_1 = crs._g1_tau_pow(N - d_bound);

        let lv = LVDigest {
            iip_x: iip_vk_x,
            iip_y: iip_vk_y,
            iip_z: iip_vk_z,
            one_idx: 3,
            mul_z_tau_2,
            d_bound,
            tau_N_minus_d_1,
        };

        MulDigest { lv, s_x, s_y, s_z }
    }
}

/// Prover for MulCircuit: given witness w = [x,y,z,1], build LV proof.
///
#[allow(non_snake_case)]
pub fn mul_prove(crs: &CRS, dg: &MulDigest, w: &MulWitness) -> MulProof {
    let w_vec = w.to_vec();

    // Three IIP proofs for selectors s_x, s_y, s_z (all over the same witness w)
    let iip_pi_x = iip_prove(crs, &dg.s_x, &w_vec);
    let iip_pi_y = iip_prove(crs, &dg.s_y, &w_vec);
    let iip_pi_z = iip_prove(crs, &dg.s_z, &w_vec);
    let nz_pi = nonzero_prove(crs, &w_vec, dg.lv.one_idx);

    // --- QAP view of the Mul demo via the “compiler” scaffold ---
    let compiled = CompiledQAP::from_mul_demo(&w_vec);
    let qap = &compiled.qap;
    let commits = QAPCommit::commit_mul(crs, qap);

    // --- MaxDeg for the IIP witness polynomial B(X) ---
    // Rebuild B(X) as interpolation of w = [x,y,z,1] on D
    let B_poly = crs.interpolate(&w_vec);
    let shift = crs.N - dg.lv.d_bound; // N - d
    let w_hat_poly = mul_by_xk(&B_poly, shift);
    let w_hat_tau_1 = crs.commit_poly_g1(w_hat_poly.coeffs());

    // Optional sanity checks
    #[cfg(debug_assertions)]
    {
        // P(1) = 0 for the QAP
        let one = Fr::from(1u32);
        let p = qap.build_p_for_mul();
        let p_at_1 = p.evaluate(&one);
        debug_assert!(p_at_1.is_zero(), "QAP check failed: P(1) != 0");

        // If x*y=z, P(X) is the zero polynomial -> [P(τ)]_1 = identity
        let gt_p =
            <Bn254 as Pairing>::pairing(commits.p_tau_1, <Bn254 as Pairing>::G2::generator());
        debug_assert!(
            gt_p.0.is_one(),
            "QAP GT check failed: [P(τ)]_1 not identity when x*y=z"
        );
    }

    let lv = LVProof {
        iip_x: iip_pi_x,
        iip_y: iip_pi_y,
        iip_z: iip_pi_z,
        nz: nz_pi,
        w: w_vec,
        p_tau_1: commits.p_tau_1,
        h_tau_1: commits.h_tau_1,
        a_tau_1: commits.a_tau_1[0],
        b_tau_1: commits.b_tau_1[0],
        c_tau_1: commits.c_tau_1[0],
        b_tau_2: commits.b_tau_2[0],
        w_hat_tau_1,
    };

    MulProof { lv }
}
