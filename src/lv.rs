// src/lv.rs
use ark_bw6_761::{BW6_761, G1Projective as G1o, G2Projective as G2o};
use ark_ec::{pairing::Pairing, Group};
use ark_groth16::{Proof as OuterProof, VerifyingKey as OuterVK};
use blake3::Hasher;

pub struct LvDigest {
    pub s_g1_basis: Vec<G1o>,
    pub s_g2_basis: Vec<G2o>,
}

pub struct LvLinearTerms {
    pub y_slots: Vec<<BW6_761 as Pairing>::TargetField>,
}

pub fn kdf_gt(gt: &<BW6_761 as Pairing>::TargetField) -> [u8; 32] {
    let mut h = Hasher::new();
    h.update(&format!("{:?}", gt).as_bytes());
    *h.finalize().as_bytes()
}

pub fn dummy_g2() -> G2o {
    G2o::generator()
}

pub fn derive_a_from_outer_proof(
    vk: &OuterVK<BW6_761>,
    proof: &OuterProof<BW6_761>,
    slots: usize,
) -> LvLinearTerms {
    let t0 = BW6_761::pairing(proof.a, proof.b).0;
    let t1 = BW6_761::pairing(vk.alpha_g1, vk.beta_g2).0;
    let t2 = BW6_761::pairing(proof.c, vk.delta_g2).0;

    let ga = if !vk.gamma_abc_g1.is_empty() {
        BW6_761::pairing(vk.gamma_abc_g1[0], vk.gamma_g2).0
    } else {
        BW6_761::pairing(vk.alpha_g1, vk.gamma_g2).0
    };

    let mut ys = vec![t0, t1, t2, ga];
    ys.resize(slots, BW6_761::pairing(vk.alpha_g1, vk.beta_g2).0);
    LvLinearTerms { y_slots: ys }
}