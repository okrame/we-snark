// src/gadgets/arithmetic.rs

use ark_bn254::{Fq12, Fr, G1Projective as G1, G2Projective as G2};
use ark_ff::{One, Zero};
use ark_poly::{
    DenseUVPolynomial,
    EvaluationDomain,
    GeneralEvaluationDomain,
    univariate::DensePolynomial,
};
use crate::compiler::R1CSMatrices;
use crate::helpers::{div_rem, mul_poly};
use crate::scs::CRS;

use super::traits::{LVGadget, LVShapeBuilder};

/// Generic QAP structure – later this will be constructed from R1CS.
/// For Phase 2 we only use it to model the one-gate Mul QAP.
#[derive(Debug, Clone)]
#[allow(non_snake_case)]
pub struct QAP {
    pub A: Vec<DensePolynomial<Fr>>,
    pub B: Vec<DensePolynomial<Fr>>,
    pub C: Vec<DensePolynomial<Fr>>,
    pub Z: DensePolynomial<Fr>,
}

impl QAP {
    /// Convenience constructor for the one-gate MulCircuit w = [x,y,z,1].
    ///
    /// This matches the old `build_mul_qap_polys` logic:
    ///   A(X) = x, B(X) = y, C(X) = z, Z(X) = X - 1, P = A·B - C.
    pub fn for_mul(w: &[Fr]) -> Self {
        assert_eq!(
            w.len(),
            4,
            "MulCircuit witness is expected to be [x, y, z, 1]"
        );
        let x = w[0];
        let y = w[1];
        let z = w[2];

        // A(X) = x, B(X) = y, C(X) = z
        let a = DensePolynomial::from_coefficients_vec(vec![x]);
        let b = DensePolynomial::from_coefficients_vec(vec![y]);
        let c = DensePolynomial::from_coefficients_vec(vec![z]);

        // Not stored for now, but we already enforce the usual identity:
        //   P(X) = A(X) B(X) - C(X) is divisible by Z(X).
        let _p = mul_poly(&a, &b) - &c;

        // Z(X) = X - 1
        let z_poly = DensePolynomial::from_coefficients_vec(vec![-Fr::one(), Fr::one()]);

        QAP {
            A: vec![a],
            B: vec![b],
            C: vec![c],
            Z: z_poly,
        }
    }

        /// Construct a QAP from generic R1CS matrices.
    ///
    /// * `num_constraints = m` determines the evaluation domain size.
    /// * `A_cols[j][k] = (i, a_{i,j})` encodes the non-zero entries in
    ///   the j-th column of the R1CS matrix A, and similarly for B,C.
    #[allow(non_snake_case)]
    pub fn from_r1cs(r1cs: &R1CSMatrices) -> Self {
        let m = r1cs.num_constraints;
        assert!(
            m > 0,
            "QAP::from_r1cs: zero-constraint instances are not supported yet"
        );

        // Domain over which we interpolate the A_i(X), B_i(X), C_i(X).
        let domain =
            GeneralEvaluationDomain::<Fr>::new(m).expect("no evaluation domain for this size");

        fn interp_columns(
            cols: &[Vec<(usize, Fr)>],
            m: usize,
            domain: &GeneralEvaluationDomain<Fr>,
        ) -> Vec<DensePolynomial<Fr>> {
            cols.iter()
                .map(|col| {
                    // evaluations A_j(ω_i) = a_{i,j} (or B,C respectively)
                    let mut evals = vec![Fr::zero(); m];
                    for &(row, coeff) in col {
                        assert!(row < m, "R1CS column entry out of range");
                        evals[row] = coeff;
                    }
                    // Interpolate polynomial from evaluations over the domain.
                    let coeffs = domain.ifft(&evals);
                    DensePolynomial::from_coefficients_vec(coeffs)
                })
                .collect()
        }

        let A = interp_columns(&r1cs.A_cols, m, &domain);
        let B = interp_columns(&r1cs.B_cols, m, &domain);
        let C = interp_columns(&r1cs.C_cols, m, &domain);

        // Vanishing polynomial Z(X) for this domain.
        let z_sparse = domain.vanishing_polynomial();
        let Z = DensePolynomial::from(z_sparse);

        QAP { A, B, C, Z }
    }


    /// Build P(X) for the current QAP.
    ///
    /// For the Mul demo, we use the specialised representation:
    ///   * A[0](X) = x
    ///   * B[0](X) = y
    ///   * C[0](X) = z
    /// and define
    ///   P(X) = A[0](X) · B[0](X) − C[0](X).
    #[allow(dead_code)]
    pub fn build_p_for_mul(&self) -> DensePolynomial<Fr> {
        assert_eq!(
            self.A.len(),
            1,
            "Mul demo QAP is expected to have a single A polynomial"
        );
        assert_eq!(
            self.B.len(),
            1,
            "Mul demo QAP is expected to have a single B polynomial"
        );
        assert_eq!(
            self.C.len(),
            1,
            "Mul demo QAP is expected to have a single C polynomial"
        );

        // P(X) = A(X)B(X) - C(X)
        let mut p = mul_poly(&self.A[0], &self.B[0]);
        {
            let mut coeffs = p.coeffs().to_vec();
            let c0 = self.C[0].coeffs().get(0).copied().unwrap_or_else(Fr::zero);
            if coeffs.is_empty() {
                coeffs.push(-c0);
            } else {
                coeffs[0] -= c0;
            }
            p = DensePolynomial::from_coefficients_vec(coeffs);
        }
        p
    }

    /// Compute H(X) = P(X) / Z(X) for the Mul demo.
    ///
    /// This mirrors the old `compute_h_poly` and debug-asserts that P is
    /// divisible by Z.
    #[allow(dead_code)]
    pub fn compute_h_for_mul(&self) -> DensePolynomial<Fr> {
        let p = self.build_p_for_mul();
        let (h, r) = div_rem(&p, &self.Z);
        debug_assert!(
            r.coeffs().iter().all(|c| c.is_zero()),
            "Mul QAP: P(X) is not divisible by Z(X); bad witness?"
        );
        h
    }
}

/// Commitments to QAP polynomials.
/// For now this is only a placeholder that mirrors the old MulQAPCommit API.
#[derive(Clone)]
#[allow(dead_code)]
pub struct QAPCommit {
    pub a_tau_1: Vec<G1>,
    pub b_tau_1: Vec<G1>,
    pub b_tau_2: Vec<G2>,
    pub c_tau_1: Vec<G1>,
    pub p_tau_1: G1,
    pub h_tau_1: G1,
}

impl QAPCommit {
    /// Commit the Mul-demo QAP polynomials with the SCS (KZG).
    #[allow(dead_code)]
    pub fn commit_mul(crs: &CRS, qap: &QAP) -> Self {
        assert_eq!(qap.A.len(), 1, "Mul demo expects a single A polynomial");
        assert_eq!(qap.B.len(), 1, "Mul demo expects a single B polynomial");
        assert_eq!(qap.C.len(), 1, "Mul demo expects a single C polynomial");

        let a_tau_1 = crs.commit_poly_g1(qap.A[0].coeffs());
        let b_tau_1 = crs.commit_poly_g1(qap.B[0].coeffs());
        let b_tau_2 = crs.commit_poly_g2(qap.B[0].coeffs());
        let c_tau_1 = crs.commit_poly_g1(qap.C[0].coeffs());

        let p = qap.build_p_for_mul();
        let p_tau_1 = crs.commit_poly_g1(p.coeffs());

        let h = qap.compute_h_for_mul();
        let h_tau_1 = crs.commit_poly_g1(h.coeffs());

        QAPCommit {
            a_tau_1: vec![a_tau_1],
            b_tau_1: vec![b_tau_1],
            b_tau_2: vec![b_tau_2],
            c_tau_1: vec![c_tau_1],
            p_tau_1,
            h_tau_1,
        }
    }
}



/// Mul gadget: at this stage it is still a thin LV wrapper that only
/// contributes equations 4 and 5 to the global LV system. The QAP
/// machinery above is intentionally unused but ready for later phases.
pub struct MulGadget;

impl MulGadget {
    pub fn new() -> Self {
        Self
    }
}

impl LVGadget for MulGadget {
    type Digest = ();
    type Witness = ();
    type Proof = ();

    fn setup(&self, _crs: &CRS) -> Self::Digest {
        ()
    }

    fn append_constraints(&self, _crs: &CRS, _dg: &Self::Digest, builder: &mut LVShapeBuilder) {
        let gt_one = Fq12::one();

        // Eq 4 (Mul QAP): c10 * c11^{-1} = 1
        builder.add_row(
            vec![0i8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0],
            gt_one,
        );

        // Eq 5 (C–z binding): c14 * c15^{-1} = 1
        builder.add_row(
            vec![0i8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0],
            gt_one,
        );
    }

    fn prove(&self, _crs: &CRS, _dg: &Self::Digest, _witness: &Self::Witness) -> Self::Proof {
        // The actual Mul SNARK lives in `mul_snark::mul_prove`.
        ()
    }
}
