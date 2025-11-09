//src/inner.rs
use ark_bls12_377::{Bls12_377, Fr};
use ark_groth16::{Groth16, Proof, ProvingKey, VerifyingKey};
use ark_relations::r1cs::ConstraintSynthesizer;
use ark_snark::SNARK;
use rand::thread_rng;

pub type InnerCurve = Bls12_377;
pub type InnerFr = Fr;

pub fn setup<C: ConstraintSynthesizer<InnerFr>>(
    circuit: C,
) -> anyhow::Result<(ProvingKey<InnerCurve>, VerifyingKey<InnerCurve>)> {
    let mut rng = thread_rng();
    let (pk, vk) = Groth16::<InnerCurve>::circuit_specific_setup(circuit, &mut rng)?;
    Ok((pk, vk))
}

pub fn prove<C: ConstraintSynthesizer<InnerFr>>(
    pk: &ProvingKey<InnerCurve>,
    circuit: C,
) -> anyhow::Result<Proof<InnerCurve>> {
    let mut rng = thread_rng();
    let prf = Groth16::<InnerCurve>::prove(pk, circuit, &mut rng)?;
    Ok(prf)
}