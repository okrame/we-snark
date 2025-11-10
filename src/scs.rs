//src/scs.rs
use ark_bn254::{Bn254, Fr, G1Projective, G2Projective};
use ark_ec::{PrimeGroup, pairing::Pairing};
use ark_ff::{Field, One, PrimeField, Zero};
use ark_poly::{
    DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain, univariate::DensePolynomial,
};
use rand::Rng;
use std::ops::Mul;

#[allow(non_snake_case)]
pub struct CRS {
    pub n: usize,                            // domain size (power of two)
    pub n_inv: Fr,                           // 1/n (y* in Construction 6 when x* = 0)
    pub g1_pows: Vec<G1Projective>,          // [tau^0]_1 .. [tau^N]_1
    pub g2_pows: Vec<G2Projective>,          // [tau^0]_2 .. [tau^N]_2
    pub N: usize,                            // max degree supported by CRS
    pub vanishing_coeffs: Vec<Fr>,           // coeffs of Z_D(X)
    pub domain: GeneralEvaluationDomain<Fr>, // D (roots of unity)
}
#[allow(non_snake_case)]
impl CRS {
    pub fn setup<R: Rng>(mut rng: R, n: usize) -> Self {
        // n must be power-of-two
        let domain = GeneralEvaluationDomain::<Fr>::new(n).expect("radix-2 domain");
        let n_inv = Fr::from(n as u64).inverse().unwrap(); // y* = 1/n at x* = 0
        let tau = Fr::from(rng.random::<u128>()); // trapdoor, local only
        // choose N >= 2n so we have indices N-n+2 and N (as in Construction 6)
        let N = 2 * n + 4;

        let mut g1_pows = Vec::with_capacity(N + 1);
        let mut g2_pows = Vec::with_capacity(N + 1);
        let g1 = <Bn254 as Pairing>::G1::generator();
        let g2 = <Bn254 as Pairing>::G2::generator();

        let mut tpow = Fr::one();
        for _ in 0..=N {
            g1_pows.push(g1.mul_bigint(tpow.into_bigint()));
            g2_pows.push(g2.mul_bigint(tpow.into_bigint()));
            tpow *= tau;
        }

        // Convert to DensePolynomial to get coefficients
        let Z_dense = DensePolynomial::from_coefficients_vec({
            let mut coeffs = vec![Fr::zero(); n + 1];
            coeffs[0] = -Fr::one();
            coeffs[n] = Fr::one();
            coeffs
        });
        let vanishing_coeffs = Z_dense.coeffs().to_vec();
        CRS {
            n,
            n_inv,
            g1_pows,
            g2_pows,
            N,
            vanishing_coeffs,
            domain,
        }
    }

    /// Commit polynomial in G1: returns [F(τ)]_1 = Σ f_j [τ^j]_1
    pub fn commit_poly_g1(&self, coeffs: &[Fr]) -> G1Projective {
        coeffs
            .iter()
            .enumerate()
            .fold(G1Projective::zero(), |acc, (j, c)| {
                if c.is_zero() {
                    acc
                } else {
                    acc + self.g1_pows[j].mul_bigint((*c).into_bigint())
                }
            })
    }

    /// Commit polynomial in G2: returns [F(τ)]_2 = Σ f_j [τ^j]_2
    pub fn commit_poly_g2(&self, coeffs: &[Fr]) -> G2Projective {
        coeffs
            .iter()
            .enumerate()
            .fold(G2Projective::zero(), |acc, (j, c)| {
                if c.is_zero() {
                    acc
                } else {
                    acc + self.g2_pows[j].mul_bigint((*c).into_bigint())
                }
            })
    }

    /// Interpolate evaluations `vals` on D to DensePolynomial coeffs
    pub fn interpolate(&self, evals: &[Fr]) -> DensePolynomial<Fr> {
        assert_eq!(evals.len(), self.n);
        // inverse FFT to get coeffs over monomial basis
        let mut v = evals.to_vec();
        self.domain.ifft_in_place(&mut v);
        DensePolynomial::from_coefficients_vec(v)
    }

    /// Multiply two polynomials (truncate/extend as needed)
    pub fn mul_poly(a: &DensePolynomial<Fr>, b: &DensePolynomial<Fr>) -> DensePolynomial<Fr> {
        a.mul(b)
    }

    pub fn div_rem(
        P: &DensePolynomial<Fr>,
        Q: &DensePolynomial<Fr>,
    ) -> (DensePolynomial<Fr>, DensePolynomial<Fr>) {
        let q = P / Q;
        let r = P - &(&q * Q);
        (q, r)
    }

    pub fn poly_from_coeffs(coeffs: Vec<Fr>) -> DensePolynomial<Fr> {
        DensePolynomial::from_coefficients_vec(coeffs)
    }

    /// Convenience: [τ^k]_2 in G2
    pub fn g2_tau_pow(&self, k: usize) -> G2Projective {
        self.g2_pows[k]
    }

    /// Convenience: [τ^k]_1 in G1
    pub fn g1_tau_pow(&self, k: usize) -> G1Projective {
        self.g1_pows[k]
    }
}
