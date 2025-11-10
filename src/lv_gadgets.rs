// src/lv_gadgets.rs
use ark_bw6_761::{BW6_761, Fr as F, G1Projective as G1, G2Projective as G2};
use ark_ec::{pairing::Pairing, Group};
use ark_ff::{Field, Zero};
use ark_std::vec::Vec;
use ark_poly::{EvaluationDomain};
use crate::scs::{Ck, Domain, xstar_ystar_fft, z_at_tau_g2};

#[derive(Clone)]
pub struct IipIndex {
    pub s_g1: Vec<G1>,
    pub aux_g1_low: Vec<Vec<G1>>,
    pub aux_g1_high: Vec<Vec<G1>>,
}

#[derive(Clone)]
pub struct IipDigest {
    pub x_star: F,
    pub y_star: F,
    pub C_g1: G1,
    pub Z_g2: G2,
}

#[derive(Clone)]
pub struct NonZeroDigest {
    pub ZS_g2: G2,
    pub s_elem: F,
}

pub struct LagrangeTable {
    pub li_coeffs: Vec<Vec<F>>,
    pub d_points: Vec<F>,
}

pub fn precompute_lagrange(dom: &Domain) -> LagrangeTable {
    let n = dom.n;
    let mut li_coeffs = vec![vec![F::ZERO; n]; n];
    for i in 0..n {
        let mut e = vec![F::ZERO; n];
        e[i] = F::ONE;
        let coeffs = dom.d.ifft(&e);
        li_coeffs[i] = coeffs;
    }
    let d_points = dom.points.clone();
    LagrangeTable { li_coeffs, d_points }
}

pub struct IIP;

impl IIP {
    pub fn auxgen(ck: &Ck, s_vec: &[F], nmax: usize) -> IipIndex {
        let n = ck.domain.n;
        let mut s_g1 = Vec::with_capacity(s_vec.len());
        let mut aux_g1_low = Vec::with_capacity(s_vec.len());
        let mut aux_g1_high = Vec::with_capacity(s_vec.len());
        
        for si in s_vec {
            s_g1.push(G1::generator() * si);
            
            let mut low = Vec::with_capacity(n);
            for j in 0..n {
                low.push(ck.g1.tau_pows[j] * si);
            }
            aux_g1_low.push(low);
            
            let mut high = Vec::with_capacity(n);
            for j in 0..n {
                let idx = nmax.saturating_sub(n) + j;
                if idx < ck.g1.tau_pows.len() {
                    high.push(ck.g1.tau_pows[idx] * si);
                } else {
                    high.push(G1::zero());
                }
            }
            aux_g1_high.push(high);
        }
        
        IipIndex { s_g1, aux_g1_low, aux_g1_high }
    }

    pub fn digest(ck: &Ck, idx: &IipIndex, lag: &LagrangeTable) -> IipDigest {
        let (x_star, y_star) = xstar_ystar_fft(&ck.domain);
        
        let c_g1 = idx.s_g1.iter().enumerate().fold(G1::zero(), |acc, (i, si_g1)| {
            let li_at_xstar: F = lag.li_coeffs[i].iter().enumerate()
                .fold(F::ZERO, |sum, (j, &coeff)| sum + coeff * x_star.pow([j as u64]));
            acc + (*si_g1 * (y_star * li_at_xstar))
        });
        
        let z_g2 = z_at_tau_g2(ck);
        
        IipDigest { x_star, y_star, C_g1: c_g1, Z_g2: z_g2 }
    }
}

pub fn qx_tau_in_exponent(
    ck: &Ck,
    idx: &IipIndex,
    w: &[F],
    lag: &LagrangeTable,
) -> G1 {
    let n = ck.domain.n;
    let (x_star, y_star) = xstar_ystar_fft(&ck.domain);
    
    let mut s_commit = G1::zero();
    for (i, si_g1) in idx.s_g1.iter().enumerate() {
        s_commit += *si_g1 * w[i];
    }
    
    let mut s_tau_low: Vec<G1> = vec![G1::zero(); n];
    for j in 0..n {
        for (k, aux_low_k) in idx.aux_g1_low.iter().enumerate() {
            if j < aux_low_k.len() {
                s_tau_low[j] += aux_low_k[j] * w[k];
            }
        }
    }
    
    let mut out = G1::zero();
    for i in 0..n {
        let denom = lag.d_points[i] - x_star;
        let inv = denom.inverse().unwrap();
        let alpha = w[i] * inv;
        let beta = -(F::ONE / y_star) * inv;
        
        for j in 0..n {
            let coeff = lag.li_coeffs[i][j] * alpha;
            if j < idx.aux_g1_low[i].len() {
                out += idx.aux_g1_low[i][j] * coeff;
            }
        }
        
        for j in 0..n {
            let coeff = lag.li_coeffs[i][j] * beta;
            out += s_tau_low[j] * coeff;
        }
    }
    out
}

pub fn qz_tau_in_exponent(
    ck: &Ck,
    idx: &IipIndex,
    w: &[F],
    lag: &LagrangeTable,
    qx_tau_g1: &G1,
) -> G1 {
    let n = ck.domain.n;
    let (x_star, y_star) = xstar_ystar_fft(&ck.domain);
    
    let b_coeffs = ck.domain.d.ifft(w);
    
    let a_tau_low = build_a_tau_shifted_low(ck, idx, lag);
    
    let mut ab_minus_sy_g1 = G1::zero();
    for j in 0..n {
        ab_minus_sy_g1 += a_tau_low[j] * b_coeffs[j];
    }
    
    let mut s_commit = G1::zero();
    for (i, si) in idx.s_g1.iter().enumerate() {
        s_commit += *si * w[i];
    }
    ab_minus_sy_g1 -= s_commit * (F::ONE / y_star);
    
    let num_tau_g1 = ab_minus_sy_g1 - (*qx_tau_g1 * x_star);
    
    deconvolve_mod_xn_minus_1(ck, idx, &num_tau_g1, lag, w)
}

fn build_a_tau_shifted_low(ck: &Ck, idx: &IipIndex, lag: &LagrangeTable) -> Vec<G1> {
    let n = ck.domain.n;
    let mut result = vec![G1::zero(); n];
    
    for j in 0..n {
        for i in 0..idx.s_g1.len() {
            for k in 0..n {
                let li_coeff = lag.li_coeffs[i][k];
                let shift_idx = (j + k) % n;
                if k < idx.aux_g1_low[i].len() {
                    result[shift_idx] += idx.aux_g1_low[i][k] * li_coeff;
                }
            }
        }
    }
    result
}

fn deconvolve_mod_xn_minus_1(
    ck: &Ck,
    idx: &IipIndex,
    num_g1: &G1,
    lag: &LagrangeTable,
    w: &[F],
) -> G1 {
    let n = ck.domain.n;
    
    let qz_coeffs_field = compute_qz_coeffs_field(ck, lag, w, n);
    
    let mut result = G1::zero();
    for (i, coeff) in qz_coeffs_field.iter().enumerate() {
        if i < idx.aux_g1_high[0].len() {
            for k in 0..idx.s_g1.len() {
                if i < idx.aux_g1_high[k].len() {
                    result += idx.aux_g1_high[k][i] * (*coeff / F::from(n as u64));
                }
            }
        }
    }
    
    result + (*num_g1 * (F::ONE / F::from(n as u64)))
}

fn compute_qz_coeffs_field(_ck: &Ck, _lag: &LagrangeTable, w: &[F], n: usize) -> Vec<F> {
    let mut coeffs = vec![F::ZERO; n];
    let inv_n = F::ONE / F::from(n as u64);
    
    for i in 0..n.min(w.len()) {
        coeffs[i] = w[i] * inv_n;
    }
    
    coeffs
}

#[derive(Clone)]
pub struct ZeroProof {
    pub q_tau_g1: G1,
}

pub fn nonzero_digest(ck: &Ck, s_elem: F) -> NonZeroDigest {
    let zs_g2 = ck.g2.tau_pows[1] - G2::generator() * s_elem;
    NonZeroDigest { ZS_g2: zs_g2, s_elem }
}

pub fn zero_prove(
    ck: &Ck,
    w: &[F],
    digest: &NonZeroDigest,
) -> ZeroProof {
    let w_at_s = w.iter().zip(ck.domain.points.iter())
        .fold(F::ZERO, |acc, (wi, di)| {
            let li_at_s = lagrange_basis_at_point(ck, di, &digest.s_elem);
            acc + (*wi * li_at_s)
        });
    
    let q_val = w_at_s;
    let q_tau_g1 = G1::generator() * q_val;
    
    ZeroProof { q_tau_g1 }
}

pub fn zero_verify(
    w_commit_g2: G2,
    proof: &ZeroProof,
    digest: &NonZeroDigest,
) -> bool {
    let lhs = BW6_761::pairing(G1::generator(), w_commit_g2).0;
    let rhs = BW6_761::pairing(proof.q_tau_g1, digest.ZS_g2).0;
    lhs == rhs
}

fn lagrange_basis_at_point(ck: &Ck, domain_point: &F, eval_point: &F) -> F {
    let mut num = F::ONE;
    let mut denom = F::ONE;
    
    for di in ck.domain.points.iter() {
        if di != domain_point {
            num *= *eval_point - di;
            denom *= *domain_point - di;
        }
    }
    
    if denom.is_zero() {
        F::ZERO
    } else {
        num / denom
    }
}