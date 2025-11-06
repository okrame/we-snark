// src/types.rs

use ark_bn254::{Fr, G1Projective, Bn254};
use ark_groth16::Proof;

/// Groth16 proof + LV "a(pi)" vector in Fr^3
pub struct LvMulProof {
    pub groth: Proof<Bn254>,
    pub lambdas: [Fr; 3],
}

/// Ciphertext format for this WE instantiation
#[allow(non_snake_case)]
pub struct WeCiphertext {
    /// instance u = z (public)
    pub u: Fr,
    /// header with S_i = g1^{s_i}
    pub S_vec: [G1Projective; 3],
    /// AES-GCM nonce
    pub nonce: [u8; 12],
    /// ciphertext bytes
    pub ct: Vec<u8>,
}
