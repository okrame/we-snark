use ark_bw6_761::Fr as OuterFr;
use ark_bls12_377::{Bls12_377, Fr as InnerFr};
use ark_groth16::{Groth16, Proof, VerifyingKey};

use ark_crypto_primitives::snark::{
    constraints::SNARKGadget
};

use ark_r1cs_std::{
    alloc::AllocVar,
    boolean::Boolean,
    eq::EqGadget,
};
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError};
use ark_relations::ns;

type InnerPairingVar = ark_r1cs_std::pairing::bls12::PairingVar<ark_bls12_377::Config>;

// Groth16 over BLS12-377
type InnerSNARK = Groth16<Bls12_377>;
// Its verification gadget
type InnerSNARKGadget = ark_groth16::constraints::Groth16VerifierGadget<Bls12_377, InnerPairingVar>;

#[derive(Clone)]
pub struct OuterVerifiesInnerCircuit {
    pub inner_vk: VerifyingKey<Bls12_377>,
    pub inner_public_inputs: Vec<InnerFr>,
    pub inner_proof: Proof<Bls12_377>,
}

impl ConstraintSynthesizer<OuterFr> for OuterVerifiesInnerCircuit {
    fn generate_constraints(self, cs: ConstraintSystemRef<OuterFr>) -> Result<(), SynthesisError> {
        // 1. Allocate inner VK as constant
        let vk_var = <InnerSNARKGadget as SNARKGadget<
            InnerFr,
            OuterFr,
            InnerSNARK,
        >>::VerifyingKeyVar::new_constant(
            ns!(cs, "vk"),
            self.inner_vk.clone(),
        )?;

        // 2. Allocate inner proof as witness
        let proof_var = <InnerSNARKGadget as SNARKGadget<
            InnerFr,
            OuterFr,
            InnerSNARK,
        >>::ProofVar::new_witness(
            ns!(cs, "proof"),
            || Ok(self.inner_proof.clone()),
        )?;

        // 3. Allocate inner public inputs as *input gadget*
        //    This wraps all the BooleanInputVar plumbing for you.
        let input_var = <InnerSNARKGadget as SNARKGadget<
            InnerFr,
            OuterFr,
            InnerSNARK,
        >>::InputVar::new_input(
            ns!(cs, "inputs"),
            || Ok(self.inner_public_inputs.clone()),
        )?;

        // 4. Call the gadget verifier
        let ok = <InnerSNARKGadget as SNARKGadget<
            InnerFr,
            OuterFr,
            InnerSNARK,
        >>::verify(&vk_var, &input_var, &proof_var)?;

        // 5. Enforce that verification result is TRUE
        ok.enforce_equal(&Boolean::constant(true))?;

        Ok(())
    }
}
