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
#[derive(Clone, Debug)]
pub enum HeaderElem { G1(G1), G2(G2) }

#[derive(Clone, Debug)]
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

fn kdf_from_gt_with_ctx(gt: &Fq12, hdr: &LVHeader, crs: &CRS, shape: &LVShape) -> [u8; 32] {
    let mut hasher = Sha256::new();
    
    // 1) GT element
    let mut gt_bytes = Vec::new();
    gt.serialize_compressed(&mut gt_bytes).unwrap();
    hasher.update(&gt_bytes);
    
    // 2) CRS context
    hasher.update(&crs.n.to_le_bytes());
    hasher.update(&crs.N.to_le_bytes());
    
    // 3) Shape matrix
    for i in 0..shape.rows {
        for j in 0..LV_NUM_COORDS {
            hasher.update(&[shape.a[i][j] as u8]);
        }
    }
    for i in 0..shape.rows {
        let mut b_bytes = Vec::new();
        shape.b[i].serialize_compressed(&mut b_bytes).unwrap();
        hasher.update(&b_bytes);
    }
    
    // 4) Header elements
    for elem in &hdr.c1 {
        let mut bytes = Vec::new();
        match elem {
            HeaderElem::G1(g) => g.serialize_compressed(&mut bytes).unwrap(),
            HeaderElem::G2(g) => g.serialize_compressed(&mut bytes).unwrap(),
        }
        hasher.update(&bytes);
    }
    
    let digest = hasher.finalize();
    let mut key = [0u8; 32];
    key.copy_from_slice(&digest);
    key
}

// binding to ct
fn compute_aad(crs: &CRS, shape: &LVShape, hdr: &LVHeader) -> Vec<u8> {
    let mut hasher = Sha256::new();
    
    hasher.update(&crs.n.to_le_bytes());
    hasher.update(&crs.N.to_le_bytes());
    
    for i in 0..shape.rows {
        for j in 0..LV_NUM_COORDS {
            hasher.update(&[shape.a[i][j] as u8]);
        }
    }
    
    for i in 0..shape.rows {
        let mut b_bytes = Vec::new();
        shape.b[i].serialize_compressed(&mut b_bytes).unwrap();
        hasher.update(&b_bytes);
    }
    
    for elem in &hdr.c1 {
        let mut bytes = Vec::new();
        match elem {
            HeaderElem::G1(g) => g.serialize_compressed(&mut bytes).unwrap(),
            HeaderElem::G2(g) => g.serialize_compressed(&mut bytes).unwrap(),
        }
        hasher.update(&bytes);
    }
    
    hasher.finalize().to_vec()
}

/// Encryptor: sample r (kept secret), compute ct1 = s·A in groups, return (header, key=H(s·b))
#[allow(non_snake_case)]
pub fn lv_make_header<R: Rng + ?Sized>(
    params: &LVPublicLinearParams,
    crs: &CRS,
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

    let hdr = LVHeader { c1 };

    // s·b in GT for KEM key (kept secret), now with context binding
    let mut B = Fq12::one();
    for i in 0..rows {
        B *= params.shape.b[i].pow(r[i].into_bigint());
    }
    let key = kdf_from_gt_with_ctx(&B, &hdr, crs, &params.shape);

    (hdr, key)
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

    Some(kdf_from_gt_with_ctx(&acc, hdr, crs, &params.shape))
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
) -> Option<Vec<u8>> {
    let key = lv_key_from_header(crs, dg, params, hdr, pi)?;
    let aad = compute_aad(crs, &params.shape, hdr);
    if aead_decrypt(key, nonce, ct, tag, &aad) {
        Some(ct.clone())
    } else {
        None
    }
}

pub fn aead_encrypt(
    crs: &CRS,
    shape: &LVShape,
    hdr: &LVHeader,
    key: [u8; 32],
    nonce_12: [u8; 12],
    plaintext: &mut Vec<u8>,
) -> Vec<u8> {
    let aad = compute_aad(crs, shape, hdr);
    let cipher = Aes256Gcm::new(&key.into());
    let nonce: &Nonce<_> = (&nonce_12).into();
    cipher
        .encrypt_in_place_detached(&nonce, &aad, plaintext)
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
    let nonce: &Nonce<_> = (&nonce_12).into();
    cipher
        .decrypt_in_place_detached(&nonce, aad, ciphertext, tag.into())
        .is_ok()
}