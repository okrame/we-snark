// src/recursion.rs
use ark_bw6_761::{BW6_761, Fr as OuterFr};
use ark_bls12_377::{Bls12_377, Fr as InnerFr};
use ark_groth16::{Groth16, Proof, ProvingKey, VerifyingKey};
use ark_relations::r1cs::ConstraintSynthesizer;
use ark_snark::SNARK;
use rand::thread_rng;

use crate::circuits::outer_verify_inner::OuterVerifiesInnerCircuit;

pub type OuterCurve = BW6_761;

pub fn setup_outer<C: ConstraintSynthesizer<OuterFr>>(
    circuit: C,
) -> anyhow::Result<(ProvingKey<OuterCurve>, VerifyingKey<OuterCurve>)> {
    let mut rng = thread_rng();
    let (pk, vk) = Groth16::<OuterCurve>::circuit_specific_setup(circuit, &mut rng)?;
    Ok((pk, vk))
}

pub fn prove_outer(
    pk_outer: &ProvingKey<OuterCurve>,
    inner_vk: VerifyingKey<Bls12_377>,
    inner_public_inputs: Vec<InnerFr>,
    inner_proof: Proof<Bls12_377>,
) -> anyhow::Result<Proof<OuterCurve>> {
    let mut rng = thread_rng();
    let circuit = OuterVerifiesInnerCircuit {
        inner_vk,
        inner_public_inputs,
        inner_proof,
    };
    let proof = Groth16::<OuterCurve>::prove(pk_outer, circuit, &mut rng)?;
    Ok(proof)
}