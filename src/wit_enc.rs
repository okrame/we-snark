//src/wit_enc.rs

/// Public LV vector b(u) in Fr^3 for the instance u = z.
///
/// In the toy construction:
///   b(u) = [1, 1, 0]
/// which does not depend on u explicitly, but conceptually b is circuit-dependent.
/// 

use ark_bn254::{Bn254, Fr, G1Projective, G2Projective};
use ark_ec::pairing::Pairing;
use ark_ec::Group;
use ark_ec::CurveGroup;
use ark_ff::{UniformRand, Zero, One};
use ark_groth16::VerifyingKey;
use aes_gcm::{Aes256Gcm, KeyInit, aead::{Aead, Payload}};
use rand::thread_rng;
use rand::RngCore;
use sha2::{Sha256, Digest};
use crate::types::WeCiphertext;

fn b_of_u(u: Fr) -> [Fr; 3] {
    [u, -u, Fr::one()]
}

/// Encryptor API:
///
/// Given:
///   - vk: verifying key (not actually used here, but realistic API)
///   - u = z: public instance
///   - msg: plaintext bytes
///
/// Returns:
///   - WeCiphertext containing (u, (S_i), nonce, ct)
///
/// The encryptor never computes a Groth16 proof.
#[allow(non_snake_case)]
pub fn we_encrypt(
    _vk: &VerifyingKey<Bn254>,
    u: Fr,
    msg: &[u8],
) -> Result<WeCiphertext, Box<dyn std::error::Error>> {
    let mut rng = thread_rng();

    // 1. Sample secret vector s in Fr^3
    let s1 = Fr::rand(&mut rng);
    let s2 = Fr::rand(&mut rng);
    let s3 = Fr::rand(&mut rng);
    let s_vec = [s1, s2, s3];
    
    fn truncate_str(s: &str, max_len: usize) -> String {
        if s.len() <= max_len * 2 {
            s.to_string()
        } else {
            format!("{}...{}", &s[..max_len], &s[s.len()-max_len..])
        }
    }
    
    println!("Encryptor samples random secret vector s:");
    println!("  s = [{}, {}, {}]", 
        truncate_str(&s1.to_string(), 8),
        truncate_str(&s2.to_string(), 8),
        truncate_str(&s3.to_string(), 8));

    // 2. Hide s in G1: S_i = g1^{s_i}
    let g1 = G1Projective::generator();
    let S_vec = [
        g1 * s1,
        g1 * s2,
        g1 * s3,
    ];
    
    println!("\nEncryptor hides this vector in G₁:");
    for (i, S_i) in S_vec.iter().enumerate() {
        let s = format!("{:?}", S_i);
        if s.len() > 60 {
            println!("  S_{} = {}...{}", i+1, &s[..40], &s[s.len()-20..]);
        } else {
            println!("  S_{} = {}", i+1, s);
        }
    }

    // 3. Compute K_exp = <s, b(u)> in Fr
    let b_vec = b_of_u(u);
    let mut K_exp = Fr::zero();
    for (s_i, b_i) in s_vec.iter().zip(b_vec.iter()) {
        K_exp += *s_i * *b_i;
    }
    
    let k_exp_str = K_exp.to_string();
    println!("\nEncryptor computes K_exp = ⟨s, b(u)⟩:");
    if k_exp_str.len() > 24 {
        println!("  K_exp = {}...{}", &k_exp_str[..12], &k_exp_str[k_exp_str.len()-12..]);
    } else {
        println!("  K_exp = {}", k_exp_str);
    }

    // 4. Turn it into a group element in G_T and derive symmetric key
    let g1K = g1 * K_exp;
    let g2 = G2Projective::generator();

    let K_gt = Bn254::pairing(g1K.into_affine(), g2.into_affine());

    // Serialize K_gt and hash to 256-bit key
    let gt_bytes = K_gt.0.to_string().into_bytes();
    let mut hasher = Sha256::new();
    hasher.update(&gt_bytes);
    let key_bytes: [u8; 32] = hasher.finalize().into();
    
    println!("\nEncryptor turns it to K ∈ G_T:");
    let k_str = format!("{:?}", K_gt);
    if k_str.len() > 60 {
        println!("  K = {}...", &k_str[..60]);
    } else {
        println!("  K = {}", k_str);
    }
    println!("\nEncryptor derives symmetric key k:");
    println!("  k = {:02x}{:02x}{:02x}{:02x}...{:02x}{:02x}{:02x}{:02x}", 
        key_bytes[0], key_bytes[1], key_bytes[2], key_bytes[3],
        key_bytes[28], key_bytes[29], key_bytes[30], key_bytes[31]);

    let cipher = Aes256Gcm::new_from_slice(&key_bytes)
        .map_err(|e| format!("AES key init error: {e}"))?;

    // 5. Generate random nonce for AES-GCM
    let mut nonce = [0u8; 12];
    rng.fill_bytes(&mut nonce);

    // 6. Encrypt
    let ct = cipher
        .encrypt(&nonce.into(), Payload { msg, aad: &[] })
        .map_err(|e| format!("Encryption failed: {e}"))?;
    
    println!("\nPublic ciphertext:");
    println!("  u = {}", u);
    println!("  |S_vec| = {} group elements", S_vec.len());
    print!("  nonce = ");
    for b in &nonce { print!("{:02x}", b); }
    println!();
    print!("  ct = ");
    for i in 0..8.min(ct.len()) { print!("{:02x}", ct[i]); }
    println!("... ({} bytes)", ct.len());

    Ok(WeCiphertext {
        u,
        S_vec,
        nonce,
        ct,
    })
}