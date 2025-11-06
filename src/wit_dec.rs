//src/wit_dec.rs
use ark_bn254::{Bn254, G1Projective, G2Projective, Fr};
use ark_ec::{pairing::Pairing, CurveGroup, Group};
use ark_ff::Zero;
use ark_snark::SNARK;
use ark_groth16::{Groth16, VerifyingKey};
use aes_gcm::{Aes256Gcm, KeyInit, aead::{Aead, Payload}};
use sha2::{Sha256, Digest};

use crate::types::{LvMulProof, WeCiphertext};
use crate::groth16_lv::derive_a_from_proof;

#[allow(non_snake_case)]
pub fn we_decrypt(
    vk: &VerifyingKey<Bn254>,
    ct: &WeCiphertext,
    lv_proof: &LvMulProof,
) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
    // 0. Build full public input vector:
    //
    //    public_inputs = [ u = z, λ1, λ2, λ3 ]
    //
    // These λ's are stored in LvMulProof and were part of the public
    // inputs when the Groth16 proof was generated.
    let public_inputs: [Fr; 4] = [
        ct.u,
        lv_proof.lambdas[0],
        lv_proof.lambdas[1],
        lv_proof.lambdas[2],
    ];
    
    println!("Decryptor receives ciphertext");


    // 1. Verify Groth16 proof against full instance (z, λ1, λ2, λ3)
    println!("\nDecryptor verifies Groth16 proof...");
    let ok = Groth16::<Bn254>::verify(vk, &public_inputs, &lv_proof.groth)?;
    if !ok {
        return Err("invalid Groth16 proof for this instance".into());
    }
    println!("   Proof is valid");

    // 2. Recover a(π) from (vk, public_inputs, proof)
    let a_vec = derive_a_from_proof(vk, &public_inputs, &lv_proof.groth)?;
    
    fn truncate_str(s: &str, max_len: usize) -> String {
        if s.len() <= max_len * 2 {
            s.to_string()
        } else {
            format!("{}...{}", &s[..max_len], &s[s.len()-max_len..])
        }
    }
    
    println!("\nDecryptor extracts a(π) from proof:");
    println!("  a(π) = [{}, {}, {}]",
        truncate_str(&a_vec[0].to_string(), 8),
        truncate_str(&a_vec[1].to_string(), 8),
        truncate_str(&a_vec[2].to_string(), 8));

    // 3. Compute T = Π S_i^{a_i} in G1
    let mut T = G1Projective::zero();
    for (S_i, a_i) in ct.S_vec.iter().zip(a_vec.iter()) {
        T += *S_i * *a_i;
    }
    
    let t_str = format!("{:?}", T);
    println!("\nDecryptor computes T = ∏ S_i^{{a_i}}:");
    if t_str.len() > 60 {
        println!("  T = {}...{}", &t_str[..40], &t_str[t_str.len()-20..]);
    } else {
        println!("  T = {}", t_str);
    }

    // 4. K' = e(T, g2)
    let g2 = G2Projective::generator();
    let K_prime = Bn254::pairing(T.into_affine(), g2.into_affine());
    
    println!("\nDecryptor computes K' = e(T, g₂):");
    let k_str = format!("{:?}", K_prime);
    if k_str.len() > 60 {
        println!("  K' = {}...", &k_str[..60]);
    } else {
        println!("  K' = {}", k_str);
    }

    // 5. Derive symmetric key from K'
    let gt_bytes = K_prime.0.to_string().into_bytes();
    let mut hasher = Sha256::new();
    hasher.update(&gt_bytes);
    let key_bytes: [u8; 32] = hasher.finalize().into();
    
    println!("\nDecryptor derives symmetric key k':");
    println!("  k' = {:02x}{:02x}{:02x}{:02x}...{:02x}{:02x}{:02x}{:02x}", 
        key_bytes[0], key_bytes[1], key_bytes[2], key_bytes[3],
        key_bytes[28], key_bytes[29], key_bytes[30], key_bytes[31]);

    let cipher = Aes256Gcm::new_from_slice(&key_bytes)
        .map_err(|e| format!("AES key init error: {e}"))?;

    // 6. Decrypt
    println!("\nDecryptor attempts AES-GCM decryption...");
    let plaintext = cipher
        .decrypt(&ct.nonce.into(), Payload { msg: &ct.ct, aad: &[] })
        .map_err(|e| format!("Decryption failed: {e}"))?;
    
    println!("  Decryption succeeded!");

    Ok(plaintext)
}