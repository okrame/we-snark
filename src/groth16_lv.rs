// src/groth16_lv.rs
//
// A concrete LV-SNARK interface on top of Groth16 for the MulCircuit:
//      relation:  x * y = z
// with z public, x,y private.
//

use ark_bn254::{Bn254, Fr};
use ark_ff::{Zero, One};
use ark_groth16::{Groth16, ProvingKey, VerifyingKey, Proof};
use ark_snark::SNARK;

use crate::circuits::simple_mul::MulCircuit;

/// LV-SNARK proof object: Groth16 proof + scalar a(pi) vector.
pub struct LvProof {
    pub groth_proof: Proof<Bn254>,
    /// a(pi) ∈ F^3 for MulCircuit
    pub a_vec: Vec<Fr>,
}

/// LV-SNARK instance object: public inputs + scalar b(u) vector.
pub struct LvInstance {
    /// public inputs for Groth16 (here: [z])
    pub public_inputs: Vec<Fr>,
    /// b(u) ∈ F^3, computable from u (= z) alone
    pub b_vec: Vec<Fr>,
}

/// A tiny helper: check that <a, b> = 0
fn inner_product_zero(a: &[Fr], b: &[Fr]) -> bool {
    assert_eq!(a.len(), b.len(), "LV vectors must have same length");
    let mut acc = Fr::zero();
    for (ai, bi) in a.iter().zip(b.iter()) {
        acc += *ai * bi;
    }
    acc.is_zero()
}

/// Construct the LV instance from vk and public input z.
pub fn lv_instance_from_vk_and_z(vk: &VerifyingKey<Bn254>, z: Fr) -> LvInstance {
    // vk is not strictly needed to compute b_vec in this toy example,
    // but we include it to make the API future-proof and consistent
    // with LV-SNARK notation (b depends on (vk, u)).
    let _ = vk; // silence unused warning

    LvInstance {
        public_inputs: vec![z],
        b_vec: vec![Fr::one(), Fr::one(), Fr::zero()],
    }
}

/// Prove in LV-SNARK form for MulCircuit:
/// - takes ProvingKey, private inputs x,y
/// - computes z = x*y
/// - generates a Groth16 proof for (x,y,z)
/// - builds a(pi) = [x*y, -z, 1]
///
/// Returns (LvProof, z).
pub fn lv_prove_mul(
    pk: &ProvingKey<Bn254>,
    x: Fr,
    y: Fr,
) -> Result<(LvProof, Fr), Box<dyn std::error::Error>> {
    let mut rng = rand::thread_rng();

    // Public output
    let z = x * y;

    // Circuit with witnesses
    let circuit = MulCircuit {
        x: Some(x),
        y: Some(y),
    };

    // Usual Groth16 proof
    let groth_proof = Groth16::<Bn254>::prove(&pk, circuit, &mut rng)?;

    // Scalar-level LV vector a(pi):
    //
    //   a(pi) = [x*y, -z, 1]
    //
    // This encodes the same R1CS relation x*y - z = 0 used by Groth16,
    // but at scalar level, suitable for Garg-style LV transformations.
    let xy = x * y;
    let a_vec = vec![
        xy,            // x*y
        -z,            // -z
        Fr::one(),     // 1
    ];

    Ok((LvProof { groth_proof, a_vec }, z))
}

/// Verify the LV-SNARK proof:
/// - standard Groth16 verification on groth_proof,
/// - plus LV check <a(pi), b(u)> = 0.
///
/// For valid proofs, both must hold.
pub fn lv_verify_mul(
    vk: &VerifyingKey<Bn254>,
    inst: &LvInstance,
    lv_proof: &LvProof,
) -> bool {
    // 1) Groth16 verification
    let groth_ok = Groth16::<Bn254>::verify(vk, &inst.public_inputs, &lv_proof.groth_proof)
        .unwrap_or(false);

    // 2) LV linear check
    let lv_ok = inner_product_zero(&lv_proof.a_vec, &inst.b_vec);

    groth_ok && lv_ok
}