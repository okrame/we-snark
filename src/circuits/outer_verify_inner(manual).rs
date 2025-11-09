use ark_bw6_761::Fr as OuterFr;
use ark_bls12_377::{Bls12_377, Fr as InnerFr};
use ark_groth16::{Proof, VerifyingKey};
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError};

use ark_r1cs_std::{
    alloc::AllocVar,
    eq::EqGadget,
    fields::nonnative::NonNativeFieldVar,
    groups::CurveVar,
    pairing::PairingVar as PairingVarTrait,
    ToBitsGadget,
};

use ark_groth16::constraints::{ProofVar, VerifyingKeyVar};

// Pairing variable for the inner curve (BLS12-377)
type InnerPairingVar = ark_r1cs_std::pairing::bls12::PairingVar<ark_bls12_377::Config>;

#[derive(Clone)]
pub struct OuterVerifiesInnerCircuit {
    pub inner_vk: VerifyingKey<Bls12_377>,
    pub inner_public_inputs: Vec<InnerFr>,
    pub inner_proof: Proof<Bls12_377>,
}

impl ConstraintSynthesizer<OuterFr> for OuterVerifiesInnerCircuit {
    fn generate_constraints(self, cs: ConstraintSystemRef<OuterFr>) -> Result<(), SynthesisError> {
        // === 1. Allocate verifying key as constant ===
        let vk_var =
            VerifyingKeyVar::<Bls12_377, InnerPairingVar>::new_constant(cs.clone(), &self.inner_vk)?;

        // === 2. Allocate proof as witness ===
        let proof_var =
            ProofVar::<Bls12_377, InnerPairingVar>::new_witness(cs.clone(), || Ok(&self.inner_proof))?;

        // === 3. Allocate inner public inputs (non-native) ===
        let mut inputs_var = Vec::new();
        for inp in &self.inner_public_inputs {
            let var =
                NonNativeFieldVar::<InnerFr, OuterFr>::new_witness(cs.clone(), || Ok(inp))?;
            inputs_var.push(var);
        }

        // === 4. Compute L = gamma_abc_g1[0] + Σ_i(x_i * gamma_abc_g1[i+1]) ===
        let mut prepared_inputs = vk_var.gamma_abc_g1[0].clone();
        for (i, input) in inputs_var.iter().enumerate() {
            let bits = input.to_bits_le()?;
            let scaled = vk_var.gamma_abc_g1[i + 1].scalar_mul_le(bits.iter())?;
            prepared_inputs += &scaled;
        }

        // === 5. Prepare group elements for pairings ===
        let prep_a = InnerPairingVar::prepare_g1(&proof_var.a)?;
        let prep_b = InnerPairingVar::prepare_g2(&proof_var.b)?;
        let prep_alpha = InnerPairingVar::prepare_g1(&vk_var.alpha_g1)?;
        let prep_beta = InnerPairingVar::prepare_g2(&vk_var.beta_g2)?;
        let prep_inputs = InnerPairingVar::prepare_g1(&prepared_inputs)?;
        let prep_gamma = InnerPairingVar::prepare_g2(&vk_var.gamma_g2)?;
        let prep_c = InnerPairingVar::prepare_g1(&proof_var.c)?;
        let prep_delta = InnerPairingVar::prepare_g2(&vk_var.delta_g2)?;

        // === 6. Compute pairing sides ===
        // LHS: e(A, B)
        let lhs = InnerPairingVar::pairing(prep_a, prep_b)?;
        // RHS: e(α, β) * e(L, γ) * e(C, δ)
        let p1 = InnerPairingVar::pairing(prep_alpha, prep_beta)?;
        let p2 = InnerPairingVar::pairing(prep_inputs, prep_gamma)?;
        let p3 = InnerPairingVar::pairing(prep_c, prep_delta)?;
        // In GT: group addition = multiplication
        let rhs = p1 + p2 + p3;

        // === 7. Enforce equality e(A,B) == e(α,β)·e(L,γ)·e(C,δ) ===
        lhs.enforce_equal(&rhs)?;

        Ok(())
    }
}
