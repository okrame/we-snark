use aes_gcm::{AeadInPlace, Aes256Gcm, KeyInit, Nonce};
use sha2::{Digest, Sha256};
use ark_ec::pairing::Pairing;
use ark_bn254::{Bn254};
use crate::iip::{IIPDigest, IIPProof};
use crate::verifier::{LVDigest, LVProof, LVShape};


/// Public parameters an encryptor will use.
pub struct LVPublicLinearParams {
    pub shape: LVShape,
    // (optionally) you can also expose the GT bases for each c_j,
    // but that's not strictly needed for the high-level A_LV·π = b_LV picture.
}

/// What the encryptor calls to obtain A_LV, b_LV.
pub fn lv_public_linear_params(crs: &CRS, dg: &LVDigest) -> LVPublicLinearParams {
    let shape = dg.linear_shape(crs);
    LVPublicLinearParams { shape }
}

/// Derive a symmetric key K = H( GT_element ) from the LV proof/digest.
/// Concretely we hash the GT element produced by the main IIP equation rearranged.
/// For simplicity we use: K_GT = (C ◦ w) / (v ◦ [y*^{-1}]_2).
pub fn derive_key_from_lv(dg: &IIPDigest, pi: &IIPProof) -> [u8;32] {
    let c_w = <Bn254 as Pairing>::pairing(dg.C, pi.w_tau_2);
    let y_inv = dg.y_star.inverse().unwrap();
    let v_scaled = pi.v_g1.mul_bigint(y_inv.into_bigint());
    let v_term = <Bn254 as Pairing>::pairing(v_scaled, <Bn254 as Pairing>::G2::generator());
    let k_gt = c_w * v_term.inverse().unwrap(); // In a valid proof: equals Q terms product in GT

    let mut bytes = Vec::new();
    k_gt.serialize_compressed(&mut bytes).unwrap();
    let key = Sha256::digest(&bytes);
    let mut out = [0u8;32];
    out.copy_from_slice(&key);
    out
}

pub fn aead_encrypt(key: [u8;32], nonce_12: [u8;12], plaintext: &mut Vec<u8>, aad: &[u8]) -> Vec<u8> {
    let cipher = Aes256Gcm::new(&key.into());
    let nonce = Nonce::from_slice(&nonce_12);
    cipher.encrypt_in_place_detached(nonce, aad, plaintext).unwrap().to_vec()
}

pub fn aead_decrypt(key: [u8;32], nonce_12: [u8;12], ciphertext: &mut Vec<u8>, tag: &[u8], aad: &[u8]) -> bool {
    let cipher = Aes256Gcm::new(&key.into());
    let nonce = Nonce::from_slice(&nonce_12);
    cipher.decrypt_in_place_detached(nonce, aad, ciphertext, tag.into()).is_ok()
}

/// Wrapper that checks LV first, then derives key, then decrypts.
pub fn decrypt_with_lv(dg: &LVDigest, pi: &LVProof, nonce: [u8;12], ct: &mut Vec<u8>, tag: &[u8], aad: &[u8]) -> Option<Vec<u8>> {
    if !crate::verifier::lv_verify(/*crs not needed here*/ &crate::scs::CRS::setup(rand::thread_rng(), dg.iip.n), dg, pi) {
        return None;
    }
    let k = derive_key_from_lv(&dg.iip, &pi.iip);
    if aead_decrypt(k, nonce, ct, tag, aad) {
        Some(ct.clone())
    } else {
        None
    }
}
