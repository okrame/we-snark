//src/we.rs
use aes_gcm::{AeadInPlace, Aes256Gcm, KeyInit, Nonce};
use sha2::{Digest, Sha256};
use ark_ff::{Field, PrimeField, Zero, One};
use ark_bn254::{Fr, Fq12, G1Projective as G1, G2Projective as G2, Bn254};
use ark_ec::pairing::Pairing;
use ark_ec::PrimeGroup;
use ark_serialize::CanonicalSerialize;
use rand::Rng;
use crate::verifier::{LVDigest, LVProof, LVShape, LV_NUM_COORDS, LVColMeta, ColSide, build_proof_side_elems};
use crate::scs::CRS;

/// LV header containing ct1 = s·A in source groups
#[derive(Clone)]
pub enum HeaderElem { G1(G1), G2(G2) }

#[derive(Clone)]
pub struct LVHeader {
    pub c1: Vec<HeaderElem>,
}

/// Public parameters an encryptor will use.
pub struct LVPublicLinearParams {
    pub shape: LVShape,
    pub cols: [LVColMeta; LV_NUM_COORDS],
}

/// What the encryptor calls to obtain A_LV, b_LV.
pub fn lv_public_linear_params(crs: &CRS, dg: &LVDigest) -> LVPublicLinearParams {
    let shape = dg.linear_shape(crs);
    let cols = dg.column_metadata(crs);
    LVPublicLinearParams { shape, cols }
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

/// Encryptor: sample r (kept secret), compute ct1 = s·A in groups, return (header, key=H(s·b))
#[allow(non_snake_case)]
pub fn lv_make_header<R: Rng + ?Sized>(
    params: &LVPublicLinearParams,
    rng: &mut R,
) -> (LVHeader, [u8; 32]) {
    let rows = params.shape.rows;

    // sample s = r (kept secret, not published)
    let mut r = Vec::with_capacity(rows);
    for _ in 0..rows {
        let mut buf = [0u8; 32];
        rng.fill(&mut buf);
        r.push(Fr::from_le_bytes_mod_order(&buf));
    }

    // α = A^T · r (field vector)
    let alpha = derive_alphas(&params.shape, &r);

    // ct1[j] = (public_base_j)^{α_j} in the appropriate source group
    let mut c1 = Vec::with_capacity(LV_NUM_COORDS);
    for j in 0..LV_NUM_COORDS {
        match params.cols[j].side {
            ColSide::ProofG1PublicG2 => {
                let base = params.cols[j].g2_pub.expect("public G2 base");
                c1.push(HeaderElem::G2(base.mul_bigint(alpha[j].into_bigint())));
            }
            ColSide::ProofG2PublicG1 => {
                let base = params.cols[j].g1_pub.expect("public G1 base");
                c1.push(HeaderElem::G1(base.mul_bigint(alpha[j].into_bigint())));
            }
        }
    }

    // s·b in GT for KEM key (kept secret)
    let mut B = Fq12::one();
    for i in 0..rows {
        B *= params.shape.b[i].pow(r[i].into_bigint());
    }
    let key = kdf_from_gt(&B);

    (LVHeader { c1 }, key)
}

/// Decryptor: derive key by pairing ct1 with proof elements to compute s·b in GT
pub fn lv_key_from_header(
    crs: &CRS,
    dg: &LVDigest,
    params: &LVPublicLinearParams,
    hdr: &LVHeader,
    pi: &LVProof,
) -> Option<[u8; 32]> {
    if hdr.c1.len() != LV_NUM_COORDS { return None; }

    if !crate::verifier::lv_verify(crs, dg, pi) { return None; }

    let proof_elems = build_proof_side_elems(crs, dg, pi)?;

    // Compute ∏_j e(proof_side_j, ct1[j]) = ∏_i b_i^{r_i} via bilinearity
    let mut acc = Fq12::one();
    for j in 0..LV_NUM_COORDS {
        match (params.cols[j].side, &hdr.c1[j], &proof_elems[j]) {
            (ColSide::ProofG1PublicG2, HeaderElem::G2(hg2), crate::verifier::ProofElem::G1(pg1)) => {
                acc *= <Bn254 as Pairing>::pairing(*pg1, *hg2).0;
            }
            (ColSide::ProofG2PublicG1, HeaderElem::G1(hg1), crate::verifier::ProofElem::G2(pg2)) => {
                acc *= <Bn254 as Pairing>::pairing(*hg1, *pg2).0;
            }
            _ => return None,
        }
    }

    Some(kdf_from_gt(&acc))
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
    let key = lv_key_from_header(crs, dg, params, hdr, pi)?;
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