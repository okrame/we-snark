// src/scs.rs
use ark_bw6_761::{Fr as F, G1Projective as G1, G2Projective as G2};
use ark_ec::Group;
use ark_ff::{Field, Zero, UniformRand};
use ark_poly::{DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain, univariate::DensePolynomial};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{vec::Vec, rand::Rng};
use anyhow::Result;

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct SrsG1 {
    pub tau_pows: Vec<G1>,
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct SrsG2 {
    pub tau_pows: Vec<G2>,
}

pub struct Domain {
    pub n: usize,
    pub d: GeneralEvaluationDomain<F>,
    pub points: Vec<F>,
}

impl Clone for Domain {
    fn clone(&self) -> Self {
        let d = GeneralEvaluationDomain::<F>::new(self.n).unwrap();
        Self {
            n: self.n,
            d,
            points: self.points.clone(),
        }
    }
}

pub struct Ck {
    pub domain: Domain,
    pub g1: SrsG1,
    pub g2: SrsG2,
}

impl Clone for Ck {
    fn clone(&self) -> Self {
        Self {
            domain: self.domain.clone(),
            g1: self.g1.clone(),
            g2: self.g2.clone(),
        }
    }
}

pub fn cgen_scs<R: Rng>(rng: &mut R, n: usize, max_deg: usize) -> Result<(F, Ck)> {
    let tau = F::rand(rng);
    let d = GeneralEvaluationDomain::<F>::new(n)
        .ok_or_else(|| anyhow::anyhow!("bad domain"))?;
    
    let mut g1 = Vec::with_capacity(max_deg + 1);
    let mut g2 = Vec::with_capacity(max_deg + 1);
    let g1_gen = G1::generator();
    let g2_gen = G2::generator();

    let mut tau_pow = F::ONE;
    for _ in 0..=max_deg {
        g1.push(g1_gen * tau_pow);
        g2.push(g2_gen * tau_pow);
        tau_pow *= tau;
    }

    let points = d.elements().collect::<Vec<_>>();
    
    Ok((tau, Ck {
        domain: Domain { n, d, points },
        g1: SrsG1 { tau_pows: g1 },
        g2: SrsG2 { tau_pows: g2 },
    }))
}

pub fn commit_g1(ck: &Ck, w: &[F]) -> G1 {
    assert_eq!(w.len(), ck.domain.n);
    interpolate_eval_in_exponent_g1(&ck.g1.tau_pows, &ck.domain, w)
}

pub fn commit_g2(ck: &Ck, w: &[F]) -> G2 {
    assert_eq!(w.len(), ck.domain.n);
    interpolate_eval_in_exponent_g2(&ck.g2.tau_pows, &ck.domain, w)
}

fn interpolate_eval_in_exponent_g1(table: &[G1], dom: &Domain, w: &[F]) -> G1 {
    let poly = DensePolynomial::from_coefficients_vec(dom.d.ifft(w));
    poly.coeffs.iter().enumerate().fold(G1::zero(), |acc, (j, c)| {
        if c.is_zero() { acc } else { acc + table[j] * c }
    })
}

fn interpolate_eval_in_exponent_g2(table: &[G2], dom: &Domain, w: &[F]) -> G2 {
    let poly = DensePolynomial::from_coefficients_vec(dom.d.ifft(w));
    poly.coeffs.iter().enumerate().fold(G2::zero(), |acc, (j, c)| {
        if c.is_zero() { acc } else { acc + table[j] * c }
    })
}

pub fn xstar_ystar_fft(dom: &Domain) -> (F, F) {
    // For multiplicative FFT domains D, Li(1) = 1/n for all i
    (F::ONE, F::ONE / F::from(dom.n as u64))
}

pub fn z_at_tau_g2(ck: &Ck) -> G2 {
    let n = ck.domain.n;
    ck.g2.tau_pows[n] - G2::generator()
}

pub fn aux_bands_for_si(ck: &Ck, si: F, n: usize, nmax: usize) -> (Vec<G1>, Vec<G1>) {
    let mut low = Vec::with_capacity(n);
    for j in 0..n {
        low.push(ck.g1.tau_pows[j] * si);
    }
    let mut high = Vec::with_capacity(n);
    for j in 0..n {
        let idx = nmax - n + j;
        high.push(ck.g1.tau_pows[idx] * si);
    }
    (low, high)
}