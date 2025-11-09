//src/we_lv.rs 
use ark_bw6_761::{BW6_761, G1Projective as G1o, Fr};
use ark_ec::{pairing::Pairing, Group};
use ark_ff::{Field, PrimeField, UniformRand};
use rand::thread_rng;

use crate::lv::{LvDigest, LvLinearTerms, kdf_gt, dummy_g2};

pub struct Ciphertext {
    pub s_g1: Vec<G1o>,
    pub nonce: [u8; 12],
    pub ct: Vec<u8>,
}

pub fn setup_lv_slots(m: usize) -> LvDigest {
    let g1 = G1o::generator();
    let g2 = dummy_g2();
    LvDigest {
        s_g1_basis: vec![g1; m],
        s_g2_basis: vec![g2; m],
    }
}

pub fn enc(ed: &LvDigest, msg: &[u8]) -> Ciphertext {
    let mut rng = thread_rng();
    let m = ed.s_g1_basis.len();
    let mut s_elems = Vec::with_capacity(m);
    for basis in ed.s_g1_basis.iter() {
        let r = Fr::rand(&mut rng);
        s_elems.push(basis.mul_bigint(r.into_bigint()));
    }
    
    let nonce_hash = blake3::hash(b"nonce");
    let mut nonce = [0u8; 12];
    nonce.copy_from_slice(&nonce_hash.as_bytes()[0..12]);
    
    Ciphertext { 
        s_g1: s_elems, 
        nonce, 
        ct: msg.to_vec() 
    }
}

pub fn dec(ct: &Ciphertext, lin: &LvLinearTerms) -> Vec<u8> {
    let g2 = dummy_g2();
    let mut acc = <BW6_761 as Pairing>::TargetField::ONE;
    
    for (s, y) in ct.s_g1.iter().zip(lin.y_slots.iter()) {
        let t = BW6_761::pairing(*s, g2);
        acc *= t.0 * y;
    }
    
    let _k = kdf_gt(&acc);
    ct.ct.clone()
}