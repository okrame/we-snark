// src/mul_snark.rs

use ark_bn254::Fr;

use crate::scs::CRS;
use crate::iip::{iip_digest, iip_prove};
use crate::nonzero::nonzero_prove;
use crate::verifier::{LVDigest, LVProof, lv_verify};

/// Fixed-size MulCircuit witness: w = [x, y, z, 1].
#[derive(Clone, Debug)]
pub struct MulWitness {
    pub x: Fr,
    pub y: Fr,
    pub z: Fr,
}

impl MulWitness {
    /// Convert to the evaluation vector [x, y, z, 1] on D.
    pub fn to_vec(&self) -> Vec<Fr> {
        vec![self.x, self.y, self.z, Fr::from(1u32)]
    }
}

/// Public parameters (vk) for the LV-SNARK.
/// For now this is just a wrapper around LVDigest + the public index s.
#[derive(Clone)]
pub struct MulDigest {
    pub lv: LVDigest,
    pub s: Vec<Fr>,
}

/// Proof object for MulCircuit: reuses LVProof as-is.
#[derive(Clone)]
pub struct MulProof {
    pub lv: LVProof,
}

impl MulDigest {
    /// Setup MulCircuit for domain size n = 4 and public index s = [0,0,1,0].
    ///
    /// - slots are [x, y, z, 1]
    /// - s selects the z-slot, so the public "output" is z
    /// - one_idx = 3 enforces w[3] = 1 via the NonZero gadget + field check
    pub fn setup(crs: &CRS) -> Self {
        assert_eq!(
            crs.n, 4,
            "MulCircuit is currently hard-coded for n=4 (slots [x,y,z,1])"
        );

        let s = vec![
            Fr::from(0u32),
            Fr::from(0u32),
            Fr::from(1u32),
            Fr::from(0u32),
        ];

        let iip_vk = iip_digest(crs, &s);
        let lv = LVDigest {
            iip: iip_vk,
            one_idx: 3,
        };

        MulDigest { lv, s }
    }
}

/// Prover for MulCircuit: given witness w = [x,y,z,1], build LV proof.
///
/// NOTE: lv_verify will re-check x*y = z and w[3] = 1 at verification time.
pub fn mul_prove(crs: &CRS, dg: &MulDigest, w: &MulWitness) -> MulProof {
    let w_vec = w.to_vec();

    let iip_pi = iip_prove(crs, &dg.s, &w_vec);
    let nz_pi  = nonzero_prove(crs, &w_vec, dg.lv.one_idx);

    let lv = LVProof {
        iip: iip_pi,
        nz:  nz_pi,
        w:   w_vec,
    };

    MulProof { lv }
}

/// Verifier wrapper of LV check + field-side mul relation.
///
/// This is the verifier for the NP relation x*y=z with w=[x,y,z,1].
pub fn mul_verify(crs: &CRS, dg: &MulDigest, pi: &MulProof) -> bool {
    lv_verify(crs, &dg.lv, &pi.lv)
}