// src/mul_snark.rs

use ark_bn254::{Bn254, Fr, G1Projective as G1, G2Projective as G2};
use ark_ec::pairing::Pairing;
use ark_ec::PrimeGroup;
use ark_ff::{One, Zero};
use ark_poly::{DenseUVPolynomial, Polynomial, univariate::DensePolynomial};

use crate::scs::CRS;
use crate::iip::{iip_digest, iip_prove};
use crate::nonzero::nonzero_prove;
use crate::verifier::{LVDigest, LVProof, lv_verify};

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
    pub s: Vec<Fr>,
}

/// Proof object for MulCircuit: reuses LVProof as-is.
#[derive(Clone)]
pub struct MulProof {
    pub lv: LVProof,
}

// ** We do not enforce any relation between these commitments and the LV gadgets yet; they’re just building blocks. **
/// QAP polynomials for the one-gate MulCircuit:
/// A(X) = x, B(X) = y, C(X) = z, Z(X) = X - 1, P(X) = A(X)B(X) - C(X).
#[derive(Clone)]
pub struct MulQAPPolys {
    pub a: DensePolynomial<Fr>,
    pub b: DensePolynomial<Fr>,
    pub c: DensePolynomial<Fr>,
    pub p: DensePolynomial<Fr>,
    pub z: DensePolynomial<Fr>,
}

/// KZG commitments to the QAP polynomials.
/// These are *not yet* integrated into the LV system; we build them under
/// debug mode to sanity-check the polynomial identities and prepare for
/// a future Mul gadget.
#[derive(Clone)]
pub struct MulQAPCommit {
    pub a_tau_1: G1,
    pub b_tau_1: G1,
    pub c_tau_1: G1,
    pub p_tau_1: G1,
    pub z_tau_2: G2,
}

/// Build QAP polynomials from the Mul witness w = [x,y,z,1].
fn build_mul_qap_polys(w: &MulWitness) -> MulQAPPolys {
    let x = w.x;
    let y = w.y;
    let z = w.z;

    // A(X) = x
    let a = DensePolynomial::from_coefficients_vec(vec![x]);

    // B(X) = y
    let b = DensePolynomial::from_coefficients_vec(vec![y]);

    // C(X) = z
    let c = DensePolynomial::from_coefficients_vec(vec![z]);

    // P(X) = A(X)B(X) - C(X)
    let mut p = CRS::mul_poly(&a, &b);
    {
        // subtract C(X) (constant z)
        let mut p_coeffs = p.coeffs().to_vec();
        if p_coeffs.is_empty() {
            p_coeffs.push(-z);
        } else {
            p_coeffs[0] -= z;
        }
        p = DensePolynomial::from_coefficients_vec(p_coeffs);
    }

    // Z(X) = X - 1
    let z_poly = DensePolynomial::from_coefficients_vec(vec![-Fr::one(), Fr::one()]);

    MulQAPPolys { a, b, c, p, z: z_poly }
}

/// Commit the QAP polynomials with the SCS (KZG).
fn commit_mul_qap(crs: &CRS, polys: &MulQAPPolys) -> MulQAPCommit {
    let a_tau_1 = crs.commit_poly_g1(polys.a.coeffs());
    let b_tau_1 = crs.commit_poly_g1(polys.b.coeffs());
    let c_tau_1 = crs.commit_poly_g1(polys.c.coeffs());
    let p_tau_1 = crs.commit_poly_g1(polys.p.coeffs());
    let z_tau_2 = crs.commit_poly_g2(polys.z.coeffs());

    MulQAPCommit {
        a_tau_1,
        b_tau_1,
        c_tau_1,
        p_tau_1,
        z_tau_2,
    }
}


impl MulDigest {
    /// Setup MulCircuit for domain size n = 4 and public index s = [0,0,1,0].
    ///
    /// - slots are [x, y, z, 1]
    /// - s selects the z-slot, so the public "output" is z
    /// - one_idx = 3 enforces w[3] = 1 via the NonZero gadget + field check
    pub fn setup(crs: &CRS) -> Self {
        assert_eq!(
            crs.n, 4,
            "MulCircuit is currently hard-coded for n=4 (slots [x,y,z,1])"
        );

        let s = vec![
            Fr::from(0u32),
            Fr::from(0u32),
            Fr::from(1u32),
            Fr::from(0u32),
        ];

        let iip_vk = iip_digest(crs, &s);
        let lv = LVDigest {
            iip: iip_vk,
            one_idx: 3,
        };

        MulDigest { lv, s }
    }
}

/// Prover for MulCircuit: given witness w = [x,y,z,1], build LV proof.
///
/// NOTE: lv_verify will re-check x*y = z and w[3] = 1 at verification time.
/// In debug builds, we also construct a tiny QAP-style view of the relation
/// and KZG-commit its polynomials as a preparation for a future Mul gadget.
pub fn mul_prove(crs: &CRS, dg: &MulDigest, w: &MulWitness) -> MulProof {
    let w_vec = w.to_vec();

    let iip_pi = iip_prove(crs, &dg.s, &w_vec);
    let nz_pi  = nonzero_prove(crs, &w_vec, dg.lv.one_idx);

    // --- QAP sanity checks (debug only) ---
    #[cfg(debug_assertions)]
    {
        let polys = build_mul_qap_polys(w);

        // Check P(1) = 0 <=> x*y - z = 0 in the field.
        let one = Fr::from(1u32);
        let p_at_1 = polys.p.evaluate(&one);
        debug_assert!(
            p_at_1.is_zero(),
            "QAP check failed: P(1) != 0, so x*y != z"
        );

        // Commit the QAP polynomials with KZG (no LV integration yet).
        let commits = commit_mul_qap(crs, &polys);

        // Optional: check at group level that P(X) is the zero polynomial
        // when x*y=z. Since P(X) is constant, P(X)=0 <=> [P(τ)]_1 is identity.
        let gt_p = <Bn254 as Pairing>::pairing(
            commits.p_tau_1,
            <Bn254 as Pairing>::G2::generator(),
        );
        debug_assert!(
            gt_p.0.is_one(),
            "QAP group check failed: [P(τ)]_1 is not the identity when x*y=z"
        );

        // We intentionally do NOT store these commitments in LVProof yet.
        // They will be used later when we extend LV_NUM_COORDS with a
        // dedicated Mul gadget.
        let _ = commits; // silence unused warning in debug
    }

    let lv = LVProof {
        iip: iip_pi,
        nz:  nz_pi,
        w:   w_vec,
    };

    MulProof { lv }
}


/// Verifier wrapper of LV check + field-side mul relation.
///
/// This is the verifier for the NP relation x*y=z with w=[x,y,z,1].
pub fn mul_verify(crs: &CRS, dg: &MulDigest, pi: &MulProof) -> bool {
    lv_verify(crs, &dg.lv, &pi.lv)
}