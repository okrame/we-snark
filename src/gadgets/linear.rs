// src/gadgets/linear.rs

use ark_bn254::{Bn254, Fq12, Fr};
use ark_ec::pairing::Pairing;
use ark_ec::PrimeGroup;
use ark_ff::One;

use crate::scs::CRS;
use crate::iip::{IIPDigest, IIPProof, iip_digest, iip_prove};
use crate::nonzero::{NonZeroProof, nonzero_prove};

use super::traits::{LVGadget, LVShapeBuilder};

// For now append_constraints is left as a stub so you don’t have to touch verifier.rs yet

/// IIP gadget: selects one linear functional <s, w>.
///
/// In your current demo you have three of these:
///   * s_x = [1,0,0,0]
///   * s_y = [0,1,0,0]
///   * s_z = [0,0,1,0]
pub struct IIPGadget {
    pub selector: Vec<Fr>,
}

impl IIPGadget {
    pub fn new(selector: Vec<Fr>) -> Self {
        Self { selector }
    }
}

impl LVGadget for IIPGadget {
    type Digest = IIPDigest;
    type Witness = Vec<Fr>;
    type Proof = IIPProof;

    fn setup(&self, crs: &CRS) -> Self::Digest {
        iip_digest(crs, &self.selector)
    }

    fn append_constraints(
        &self,
        _crs: &CRS,
        _dg: &Self::Digest,
        builder: &mut LVShapeBuilder,
    ) {
        // These three rows are exactly Eq 0–2 from LVDigest::linear_shape.
        let one = Fq12::one();

        // Eq 0: c0 * c1^{-1} * c2^{-1} * c3^{-1} = 1
        builder.add_row(
            vec![
                1, -1, -1, -1, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0,
            ],
            one,
        );

        // Eq 1: c4 * c5^{-1} = 1
        builder.add_row(
            vec![
                0, 0, 0, 0, 1, -1, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0,
            ],
            one,
        );

        // Eq 2: c6 * c7^{-1} = 1
        builder.add_row(
            vec![
                0, 0, 0, 0, 0, 0, 1, -1, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0,
            ],
            one,
        );
    }

    fn prove(
        &self,
        crs: &CRS,
        _dg: &Self::Digest,
        witness: &Self::Witness,
    ) -> Self::Proof {
        iip_prove(crs, &self.selector, witness)
    }
}


/// NonZero gadget: enforces that w[idx_one] == 1.
pub struct NonZeroGadget {
    pub idx_one: usize,
}

impl NonZeroGadget {
    pub fn new(idx_one: usize) -> Self {
        Self { idx_one }
    }
}

impl LVGadget for NonZeroGadget {
    type Digest = ();
    type Witness = Vec<Fr>;
    type Proof = NonZeroProof;

    fn setup(&self, _crs: &CRS) -> Self::Digest {
        ()
    }

    fn append_constraints(
        &self,
        _crs: &CRS,
        _dg: &Self::Digest,
        builder: &mut LVShapeBuilder,
    ) {
        // Eq 3: c8 * c9^{-1} = e(g1,g2)
        let g1 = <Bn254 as Pairing>::G1::generator();
        let g2 = <Bn254 as Pairing>::G2::generator();
        let gt_const: Fq12 = <Bn254 as Pairing>::pairing(g1, g2).0;

        builder.add_row(
            vec![
                0, 0, 0, 0, 0, 0, 0, 0, 1, -1,
                0, 0, 0, 0, 0, 0, 0, 0,
            ],
            gt_const,
        );
    }

    fn prove(
        &self,
        crs: &CRS,
        _dg: &Self::Digest,
        witness: &Self::Witness,
    ) -> Self::Proof {
        nonzero_prove(crs, witness, self.idx_one)
    }
}


/// Placeholder for the MaxDeg gadget that will enforce deg ≤ d_bound.
///
/// Right now the logic lives implicitly in:
///   * `LVDigest { d_bound, tau_N_minus_d_1 }`
///   * `MulProof::w_hat_tau_1`
///   * LV equation “c16 * c17^{-1} = 1”.
pub struct MaxDegGadget {
    pub d_bound: usize,
}

impl MaxDegGadget {
    pub fn new(d_bound: usize) -> Self {
        Self { d_bound }
    }
}

impl LVGadget for MaxDegGadget {
    type Digest = ();
    type Witness = Vec<Fr>;
    type Proof = ();

    fn setup(&self, _crs: &CRS) -> Self::Digest {
        ()
    }

    fn append_constraints(
        &self,
        _crs: &CRS,
        _dg: &Self::Digest,
        builder: &mut LVShapeBuilder,
    ) {
        // Eq 6 (MaxDeg for B): c16 * c17^{-1} = 1
        let one = Fq12::one();

        builder.add_row(
            vec![
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 1, -1,
            ],
            one,
        );
    }

    fn prove(
        &self,
        _crs: &CRS,
        _dg: &Self::Digest,
        _witness: &Self::Witness,
    ) -> Self::Proof {
        // Mul-related MaxDeg proof stays in mul_snark for now; this is LV-only.
        ()
    }
}

