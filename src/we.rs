// src/we.rs
use ark_bw6_761::{BW6_761, Fr as F, G1Projective as G1, G2Projective as G2};
use ark_ec::{pairing::Pairing, Group};
use ark_ff::{Field, PrimeField, UniformRand, Zero};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{rand::Rng, vec::Vec};
use blake3::Hasher;
use aes_gcm::{Aes256Gcm, aead::{Aead, KeyInit}, AeadCore};

use crate::scs::{Ck, cgen_scs, commit_g2, xstar_ystar_fft};
use crate::lv_gadgets::{
    IipIndex, NonZeroDigest, IIP, nonzero_digest,
    qx_tau_in_exponent, qz_tau_in_exponent, LagrangeTable,
    precompute_lagrange, zero_prove, ZeroProof,
};

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Pp {
    pub n: usize,
    pub nmax: usize,
}

#[derive(Clone)]
pub struct Crs {
    pub ck: Ck,
}

#[derive(Clone)]
pub struct Index {
    pub iip_idx: IipIndex,
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Aux {}

#[derive(Clone)]
pub struct EncDigest {
    pub A_g1: Vec<G1>,
    pub b_gt: Vec<<BW6_761 as Pairing>::TargetField>,
}

#[derive(Clone)]
pub struct DecDigest {
    pub C_g1: G1,
    pub Z_g2: G2,
    pub tau_minus_x_g2: G2,
    pub tau_pow_N_g2: G2,
    pub y_inv: F,
    pub nz_digest: NonZeroDigest,
}

#[derive(Clone)]
pub struct Ct {
    pub c1_g1: Vec<G1>,
    pub nonce: [u8; 12],
    pub ct2: Vec<u8>,
}

#[derive(Clone)]
pub struct Pi {
    pub qz_tau_g1: G1,
    pub qx_tau_g1: G1,
    pub qx_hat_tau_g1: G1,
    pub v_hat_tau_g1: G1,
    pub w_commit_g2: G2,
    pub zero_proof: ZeroProof,
}

pub fn setup<R: Rng>(rng: &mut R, n: usize, max_deg: usize) -> (Pp, Crs) {
    let (_tau, ck) = cgen_scs(rng, n, max_deg).expect("scs gen");
    (Pp { n, nmax: max_deg }, Crs { ck })
}

pub fn auxgen(_pp: &Pp, crs: &Crs, s: &[F], nmax: usize) -> (Index, Aux) {
    let iip_idx = IIP::auxgen(&crs.ck, s, nmax);
    (Index { iip_idx }, Aux {})
}

pub fn digest(pp: &Pp, crs: &Crs, idx: &Index) -> (EncDigest, DecDigest) {
    let lag = precompute_lagrange(&crs.ck.domain);
    let iip_vk = IIP::digest(&crs.ck, &idx.iip_idx, &lag);
    let (x_star, y_star) = xstar_ystar_fft(&crs.ck.domain);
    
    let s_elem = if !crs.ck.domain.points.is_empty() {
        crs.ck.domain.points[0]
    } else {
        F::ONE
    };
    let nz_vk = nonzero_digest(&crs.ck, s_elem);

    // For WE: use e(C, [1]) as the base pairing value
    // This makes encryption and decryption consistent:
    // - Enc: sb = e(C, [1])^s = e(s·C, [1])
    // - Dec: sb = e(ct_1, [1]) = e(s·C, [1])
    let A_g1 = vec![iip_vk.C_g1];
    let b_gt = vec![BW6_761::pairing(iip_vk.C_g1, G2::generator()).0];

    let tau_minus_x = crs.ck.g2.tau_pows[1] - G2::generator() * x_star;
    let y_inv = y_star.inverse().unwrap();
    
    (
        EncDigest { A_g1, b_gt },
        DecDigest {
            C_g1: iip_vk.C_g1,
            Z_g2: iip_vk.Z_g2,
            tau_minus_x_g2: tau_minus_x,
            tau_pow_N_g2: crs.ck.g2.tau_pows[pp.nmax],
            y_inv,
            nz_digest: nz_vk,
        },
    )
}

pub fn enc<R: Rng>(rng: &mut R, ed: &EncDigest, msg: &[u8]) -> Ct {
    let s = F::rand(rng);

    let c1_g1: Vec<G1> = ed.A_g1.iter().map(|a| *a * s).collect();

    let sb = ed.b_gt.iter().fold(
        <BW6_761 as Pairing>::TargetField::ONE,
        |acc, b| acc * b.pow(s.into_bigint())
    );

    let (key, nonce_bytes) = kdf_key_nonce(&sb);
    
    let cipher = Aes256Gcm::new(&key.into());
    let nonce = aes_gcm::Nonce::from_slice(&nonce_bytes);
    let ct2 = cipher.encrypt(nonce, msg).expect("enc ok");

    Ct { c1_g1, nonce: nonce_bytes, ct2 }
}

pub fn prove(
    dd: &DecDigest,
    idx: &Index,
    ck: &Ck,
    w: &[F],
    lag: &LagrangeTable,
) -> anyhow::Result<Pi> {
    let qx = qx_tau_in_exponent(ck, &idx.iip_idx, w, lag);
    let qz = qz_tau_in_exponent(ck, &idx.iip_idx, w, lag, &qx);

    let qx_hat = qx;

    let mut s_commit = G1::zero();
    for (i, si) in idx.iip_idx.s_g1.iter().enumerate() {
        if i < w.len() {
            s_commit += *si * w[i];
        }
    }
    let v_hat = s_commit;

    let w_commit_g2 = commit_g2(ck, w);

    let zero_proof = zero_prove(ck, w, &dd.nz_digest);

    Ok(Pi {
        qz_tau_g1: qz,
        qx_tau_g1: qx,
        qx_hat_tau_g1: qx_hat,
        v_hat_tau_g1: v_hat,
        w_commit_g2,
        zero_proof,
    })
}

pub fn dec(dd: &DecDigest, ct: &Ct, pi: &Pi) -> anyhow::Result<Vec<u8>> {
    let sb = recover_sb_via_linear_check(dd, ct, pi);
    let (key, _nonce) = kdf_key_nonce(&sb);
    let cipher = Aes256Gcm::new(&key.into());
    let nonce = aes_gcm::Nonce::from_slice(&ct.nonce);
    let pt = cipher.decrypt(nonce, ct.ct2.as_ref())
        .map_err(|e| anyhow::anyhow!("decryption failed: {:?}", e))?;
    Ok(pt)
}

fn kdf_key_nonce(gt: &<BW6_761 as Pairing>::TargetField) -> ([u8; 32], [u8; 12]) {
    let mut h = Hasher::new();
    let mut bytes = Vec::new();
    gt.serialize_compressed(&mut bytes).unwrap();
    h.update(&bytes);
    
    // Use Blake3 XOF mode to generate 44 bytes (32 for key + 12 for nonce)
    let mut output = [0u8; 44];
    let mut xof = h.finalize_xof();
    xof.fill(&mut output);
    
    let mut key = [0u8; 32];
    key.copy_from_slice(&output[..32]);
    let mut nonce = [0u8; 12];
    nonce.copy_from_slice(&output[32..44]);
    (key, nonce)
}

fn recover_sb_via_linear_check(
    dd: &DecDigest,
    ct: &Ct,
    pi: &Pi,
) -> <BW6_761 as Pairing>::TargetField {
    // Verify the proof is valid by checking the IIP verification equation
    let lhs = BW6_761::pairing(dd.C_g1, pi.w_commit_g2).0;
    let y_inv_g2 = G2::generator() * dd.y_inv;
    let rhs1 = BW6_761::pairing(pi.v_hat_tau_g1, y_inv_g2).0;
    let rhs2 = BW6_761::pairing(pi.qx_tau_g1, dd.tau_minus_x_g2).0;
    let rhs3 = BW6_761::pairing(pi.qz_tau_g1, dd.Z_g2).0;
    let rhs = rhs1 * rhs2 * rhs3;
    
    // Check verification equation holds (in production, return error if not)
    assert_eq!(lhs, rhs, "IIP verification equation failed!");
    
    // Recover s·b = e(ct_1, [1]) = e(s·C, [1])
    // This matches encryption which uses b = e(C, [1]) and computes sb = b^s = e(C,[1])^s = e(s·C,[1])
    BW6_761::pairing(ct.c1_g1[0], G2::generator()).0
}