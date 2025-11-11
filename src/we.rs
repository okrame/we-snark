//src/we.rs
use aes_gcm::{AeadInPlace, Aes256Gcm, KeyInit, Nonce};
use sha2::{Digest, Sha256};
use ark_ff::{Field, PrimeField, Zero, One};
use ark_bn254::{Fr, Fq12};
use ark_serialize::CanonicalSerialize;
use rand::Rng;
use crate::verifier::{LVDigest, LVProof, LVShape, LV_NUM_COORDS};
use crate::scs::CRS;

/// LV header containing random linear combination coefficients
#[derive(Clone)]
pub struct LVHeader {
    pub r: Vec<Fr>,
}

/// Public parameters an encryptor will use.
pub struct LVPublicLinearParams {
    pub shape: LVShape,
}

/// What the encryptor calls to obtain A_LV, b_LV.
pub fn lv_public_linear_params(crs: &CRS, dg: &LVDigest) -> LVPublicLinearParams {
    let shape = dg.linear_shape(crs);
    LVPublicLinearParams { shape }
}

fn gt_pow(base: &Fq12, exp: &Fr) -> Fq12 {
    base.pow(exp.into_bigint())
}

fn derive_alphas(shape: &LVShape, r: &[Fr]) -> [Fr; LV_NUM_COORDS] {
    let mut alpha = [Fr::zero(); LV_NUM_COORDS];
    for i in 0..shape.rows {
        let ri = r[i];
        for j in 0..LV_NUM_COORDS {
            match shape.a[i][j] {
                1  => { alpha[j] += ri; }
                -1 => { alpha[j] -= ri; }
                _  => {}
            }
        }
    }
    alpha
}

fn kdf_from_gt(gt: &Fq12) -> [u8; 32] {
    let mut bytes = Vec::new();
    gt.serialize_compressed(&mut bytes).unwrap();
    let digest = Sha256::digest(&bytes);
    let mut out = [0u8; 32];
    out.copy_from_slice(&digest);
    out
}

/// Encryptor: sample r, compute B = Π b_i^{r_i}, return (header, key=H(B))
#[allow(non_snake_case)]
pub fn lv_make_header<R: Rng + ?Sized>(
    params: &LVPublicLinearParams,
    rng: &mut R,
) -> (LVHeader, [u8; 32]) {
    let rows = params.shape.rows;
    let mut r = Vec::with_capacity(rows);

    for _ in 0..rows {
        let mut buf = [0u8; 32];
        rng.fill(&mut buf);
        let ri = Fr::from_le_bytes_mod_order(&buf);
        r.push(ri);
    }

    // RHS: B = Π_i b_i^{r_i}
    let mut B = Fq12::one();
    for i in 0..rows {
        B *= gt_pow(&params.shape.b[i], &r[i]);
    }

    let key = kdf_from_gt(&B);
    (LVHeader { r }, key)
}

/// Decryptor: recompute coords from π, aggregate with α = r·A, check vs public RHS, derive key.
pub fn lv_key_from_header(
    crs: &CRS,
    dg: &LVDigest,
    shape: &LVShape,
    hdr: &LVHeader,
    pi: &LVProof,
) -> Option<[u8; 32]> {
    if hdr.r.len() != shape.rows { return None; }

    if !crate::verifier::lv_verify(crs, dg, pi) { return None; }

    let coords = crate::verifier::build_lv_coords(crs, dg, pi)?;

    let alpha = derive_alphas(shape, &hdr.r);

    let mut lhs = Fq12::one();
    for j in 0..LV_NUM_COORDS {
        lhs *= gt_pow(&coords.0[j], &alpha[j]);
    }

    let mut rhs = Fq12::one();
    for i in 0..shape.rows {
        rhs *= gt_pow(&shape.b[i], &hdr.r[i]);
    }

    if lhs != rhs { return None; }
    Some(kdf_from_gt(&rhs))
}

pub fn decrypt_with_lv_header(
    crs: &CRS,
    dg: &LVDigest,
    params: &LVPublicLinearParams,
    hdr: &LVHeader,
    pi: &LVProof,
    nonce: [u8; 12],
    ct: &mut Vec<u8>,
    tag: &[u8],
    aad: &[u8],
) -> Option<Vec<u8>> {
    let key = lv_key_from_header(crs, dg, &params.shape, hdr, pi)?;
    if aead_decrypt(key, nonce, ct, tag, aad) {
        Some(ct.clone())
    } else {
        None
    }
}

pub fn aead_encrypt(
    key: [u8; 32],
    nonce_12: [u8; 12],
    plaintext: &mut Vec<u8>,
    aad: &[u8],
) -> Vec<u8> {
    let cipher = Aes256Gcm::new(&key.into());
    let nonce = Nonce::clone_from_slice(&nonce_12);
    cipher
        .encrypt_in_place_detached(&nonce, aad, plaintext)
        .unwrap()
        .to_vec()
}

pub fn aead_decrypt(
    key: [u8; 32],
    nonce_12: [u8; 12],
    ciphertext: &mut Vec<u8>,
    tag: &[u8],
    aad: &[u8],
) -> bool {
    let cipher = Aes256Gcm::new(&key.into());
    let nonce = Nonce::clone_from_slice(&nonce_12);
    cipher
        .decrypt_in_place_detached(&nonce, aad, ciphertext, tag.into())
        .is_ok()
}