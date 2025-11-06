// src/groth16_lv.rs

use ark_bn254::{Bn254, Fr};
use ark_groth16::{Groth16, ProvingKey, Proof, VerifyingKey};
use ark_snark::SNARK;
use rand::thread_rng; 
use ark_ff::One;

use crate::circuits::simple_mul::MulCircuit;
use crate::types::LvMulProof;


impl LvMulProof {
    // helper: expose number of coordinates
    //pub const dim: usize = 3;
}

/// Prover side: generates a Groth16 proof and the LV coordinates (λ1, λ2, λ3)
///
/// - pk: Groth16 proving key for MulCircuit
/// - x, y: private inputs (witness)
///
/// Returns: (LvMulProof, z)
///   where z = x * y is the public output / instance u.
pub fn lv_prove_mul(
    pk: &ProvingKey<Bn254>,
    x: Fr,
    y: Fr,
) -> Result<(LvMulProof, Fr), Box<dyn std::error::Error>> {
    let mut rng = thread_rng();

    // Public output z = x * y
    let z = x * y;

    // Define λ's as functions of (x, y, z):
    //   λ1 = x*y = z
    //   λ2 = -z
    //   λ3 = 1
    let lambda1 = z;
    let lambda2 = -z;
    let lambda3 = Fr::one();

    // Circuit with witnesses and public inputs filled in
    let circuit = MulCircuit {
        x: Some(x),
        y: Some(y),
        z: Some(z),
        lambda1: Some(lambda1),
        lambda2: Some(lambda2),
        lambda3: Some(lambda3),
    };

    // Standard Groth16 proof
    let groth_proof = Groth16::<Bn254>::prove(pk, circuit, &mut rng)?;

    Ok((
        LvMulProof {
            groth: groth_proof,
            lambdas: [lambda1, lambda2, lambda3],
        },
        z,
    ))
}


/// Linearly verifiable map a(π) = f(vk, u, π) ∈ Fr^3 for MulCircuit.
///
/// In this instantiation, we have extended the public input of the circuit to:
///   public_inputs = [ z, λ1, λ2, λ3 ]
///
/// The circuit enforces:
///   z = x*y
///   λ1 = x*y
///   λ2 = -z
///   λ3 = 1
///
/// So we define:
///   a(π) = [λ1, λ2, λ3] = [x*y, -z, 1].
///
/// The dependence on π is via the Groth16 proof: only a valid proof can make
/// these public inputs consistent with some witness (x,y).
pub fn derive_a_from_proof(
    _vk: &VerifyingKey<Bn254>,
    public_inputs: &[Fr],    
    _proof: &Proof<Bn254>,
) -> Result<[Fr; 3], Box<dyn std::error::Error>> {
    if public_inputs.len() != 4 {
        return Err("derive_a_from_proof: expected 4 public inputs [z, λ1, λ2, λ3]".into());
    }

    let lambda1 = public_inputs[1];
    let lambda2 = public_inputs[2];
    let lambda3 = public_inputs[3];

    Ok([lambda1, lambda2, lambda3])
}