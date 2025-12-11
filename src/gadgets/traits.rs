// src/gadgets/traits.rs

use ark_bn254::Fq12;
use crate::scs::CRS;

/// Helper for building the global LV linear system A_LV · π = b_LV
/// in a gadget-agnostic way.
#[derive(Default)]
pub struct LVShapeBuilder {
    /// Coefficients a[i][j] ∈ {-1, 0, 1} for equation i, coordinate j.
    pub a: Vec<Vec<i8>>,
    /// Right-hand sides b[i] ∈ GT for equation i.
    pub b: Vec<Fq12>,
}

impl LVShapeBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    /// Append one equation: ∏_j c_j^{coeffs[j]} = rhs.
    pub fn add_row(&mut self, coeffs: Vec<i8>, rhs: Fq12) {
        if !self.a.is_empty() {
            assert_eq!(
                self.a[0].len(),
                coeffs.len(),
                "all LV rows must have the same number of columns"
            );
        }
        self.a.push(coeffs);
        self.b.push(rhs);
    }

    pub fn rows(&self) -> usize {
        self.a.len()
    }

    pub fn cols(&self) -> usize {
        self.a.first().map(|r| r.len()).unwrap_or(0)
    }
}

/// Abstraction over a “gadget” in the Garg LV framework.
///
/// Each gadget:
///   * has a public digest (vk),
///   * knows how to append its LV rows to the global system,
///   * knows how to produce a gadget-local proof from a witness.
pub trait LVGadget {
    type Digest;
    type Witness;
    type Proof;

    /// Build the public digest (vk) for this gadget from the CRS.
    fn setup(&self, crs: &CRS) -> Self::Digest;

    /// Append this gadget's LV equations to the global system.
    ///
    /// For this to work we need `LVShape` to be made dynamic.
    fn append_constraints(
        &self,
        crs: &CRS,
        dg: &Self::Digest,
        builder: &mut LVShapeBuilder,
    );

    /// Build a gadget-local proof from a concrete witness.
    fn prove(
        &self,
        crs: &CRS,
        dg: &Self::Digest,
        witness: &Self::Witness,
    ) -> Self::Proof;
}
