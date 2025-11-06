// src/circuits/simple_mul.rs
use ark_ff::{PrimeField};
use ark_relations::lc;
use ark_relations::r1cs::{
    ConstraintSynthesizer, ConstraintSystemRef, SynthesisError, Variable,
};

#[derive(Clone)]
pub struct MulCircuit<F: PrimeField> {
    pub x: Option<F>,
    pub y: Option<F>,
    pub z: Option<F>,
    pub lambda1: Option<F>,
    pub lambda2: Option<F>,
    pub lambda3: Option<F>,
}

impl<F: PrimeField> ConstraintSynthesizer<F> for MulCircuit<F> {
    fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError> {
        // Witness variables x, y
        let x = cs.new_witness_variable(|| self.x.ok_or(SynthesisError::AssignmentMissing))?;

        let y = cs.new_witness_variable(|| self.y.ok_or(SynthesisError::AssignmentMissing))?;

        // Public inputs: z, lambda1, lambda2, lambda3
        let z = cs.new_input_variable(|| self.z.ok_or(SynthesisError::AssignmentMissing))?;

        let lambda1 =
            cs.new_input_variable(|| self.lambda1.ok_or(SynthesisError::AssignmentMissing))?;

        let lambda2 =
            cs.new_input_variable(|| self.lambda2.ok_or(SynthesisError::AssignmentMissing))?;

        let lambda3 =
            cs.new_input_variable(|| self.lambda3.ok_or(SynthesisError::AssignmentMissing))?;

        let one: Variable = Variable::One;

        // 1. Enforce z = x * y
        //
        //    (x) * (y) = z
        cs.enforce_constraint(lc!() + x, lc!() + y, lc!() + z)?;

        // 2. Enforce lambda1 = z
        //
        //    (lambda1) * 1 = z
        cs.enforce_constraint(
            lc!() + lambda1,
            lc!() + (F::one(), one),
            lc!() + z,
        )?;

        // 3. Enforce lambda2 = -z
        //
        //    (lambda2) * 1 = -z
        cs.enforce_constraint(
            lc!() + lambda2,
            lc!() + (F::one(), one),
            lc!() - z,
        )?;

        // 4. Enforce lambda3 = 1
        //
        //    (lambda3) * 1 = 1
        cs.enforce_constraint(
            lc!() + lambda3,
            lc!() + (F::one(), one),
            lc!() + (F::one(), one),
        )?;

        Ok(())
    }
}