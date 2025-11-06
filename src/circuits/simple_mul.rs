// src/circuits/simple_mul.rs
use ark_ff::PrimeField;
use ark_relations::lc;
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystemRef, SynthesisError};

#[derive(Clone)]
pub struct MulCircuit<F: PrimeField> {
    pub x: Option<F>,
    pub y: Option<F>,
}

impl<F: PrimeField> ConstraintSynthesizer<F> for MulCircuit<F> {
    fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError> {
        // Allocate private witnesses for x and y
        let x_var = cs.new_witness_variable(|| self.x.ok_or(SynthesisError::AssignmentMissing))?;
        let y_var = cs.new_witness_variable(|| self.y.ok_or(SynthesisError::AssignmentMissing))?;

        // Compute z = x * y (for witness assignment)
        let z_value = self.x.and_then(|x| self.y.map(|y| x * y));

        // Allocate public input for z
        let z_var = cs.new_input_variable(|| z_value.ok_or(SynthesisError::AssignmentMissing))?;

        // Enforce the constraint: x * y == z
        cs.enforce_constraint(lc!() + x_var, lc!() + y_var, lc!() + z_var)?;

        Ok(())
    }
}