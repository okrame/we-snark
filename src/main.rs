mod scs;
mod iip;
mod nonzero;
mod verifier;
mod we;
mod mul_snark;
mod helpers;

use ark_bn254::Fr;
use ark_serialize::CanonicalSerialize;
use rand::{rng, Rng};
use std::time::Instant;

use scs::CRS;
use we::{aead_encrypt, decrypt_with_lv_header};
use mul_snark::{MulDigest, MulWitness, mul_prove};
use crate::verifier::{lv_verify};


fn serialized_size<T: CanonicalSerialize>(t: &T) -> usize {
    t.serialized_size(ark_serialize::Compress::No)
}


fn main() {
    let mut rng = rng();

    // --- Parameters ---
    // Domain size n = 4: slots [x, y, z, 1]
    let n = 4;
    let crs = CRS::setup(&mut rng, n);

    let x = Fr::from(12u32);
    let y = Fr::from(17u32);
    let z = x * y;

    let w = MulWitness { x, y, z };

    let dg = MulDigest::setup(&crs, z);
    let pi = mul_prove(&crs, &dg, &w);

    // sanity check
    assert!(lv_verify(&crs, &dg.lv, &pi.lv));

    println!("\n=== SIZE MEASUREMENTS (bytes) ===");
    
    // Public parameters
    let crs_size = serialized_size(&crs.g1_pows) + serialized_size(&crs.g2_pows);
    println!("CRS (g1_pows + g2_pows): {}", crs_size);
    
    // Digest size: manually calculate from components
    let digest_size = serialized_size(&dg.lv.iip_x.C) + serialized_size(&dg.lv.iip_x.Z_tau_2) + 
                      serialized_size(&dg.lv.iip_x.tau_2) + serialized_size(&dg.lv.iip_x.tau_N_minus_n_plus_2_2) + 
                      serialized_size(&dg.lv.iip_x.tau_N_2) +
                      serialized_size(&dg.lv.iip_y.C) + serialized_size(&dg.lv.iip_y.Z_tau_2) + 
                      serialized_size(&dg.lv.iip_y.tau_2) + serialized_size(&dg.lv.iip_y.tau_N_minus_n_plus_2_2) + 
                      serialized_size(&dg.lv.iip_y.tau_N_2) +
                      serialized_size(&dg.lv.iip_z.C) + serialized_size(&dg.lv.iip_z.Z_tau_2) + 
                      serialized_size(&dg.lv.iip_z.tau_2) + serialized_size(&dg.lv.iip_z.tau_N_minus_n_plus_2_2) + 
                      serialized_size(&dg.lv.iip_z.tau_N_2) +
                      serialized_size(&dg.lv.mul_z_tau_2) + serialized_size(&dg.lv.instance_z) + 
                      serialized_size(&dg.lv.tau_N_minus_d_1);
    println!("Digest (verification key): {}", digest_size);
    
    // Witness
    let witness_size = serialized_size(&w.x) + serialized_size(&w.y) + serialized_size(&w.z);
    println!("Witness (x, y, z): {}", witness_size);
    
    // Proof size: manually calculate from components
    let proof_size = serialized_size(&pi.lv.iip_x.w_tau_2) + serialized_size(&pi.lv.iip_x.v_g1) +
                     serialized_size(&pi.lv.iip_x.QZ_tau_1) + serialized_size(&pi.lv.iip_x.QX_tau_1) +
                     serialized_size(&pi.lv.iip_x.QX_hat_tau_1) + serialized_size(&pi.lv.iip_x.v_hat_tau_1) +
                     serialized_size(&pi.lv.iip_y.w_tau_2) + serialized_size(&pi.lv.iip_y.v_g1) +
                     serialized_size(&pi.lv.iip_y.QZ_tau_1) + serialized_size(&pi.lv.iip_y.QX_tau_1) +
                     serialized_size(&pi.lv.iip_y.QX_hat_tau_1) + serialized_size(&pi.lv.iip_y.v_hat_tau_1) +
                     serialized_size(&pi.lv.iip_z.w_tau_2) + serialized_size(&pi.lv.iip_z.v_g1) +
                     serialized_size(&pi.lv.iip_z.QZ_tau_1) + serialized_size(&pi.lv.iip_z.QX_tau_1) +
                     serialized_size(&pi.lv.iip_z.QX_hat_tau_1) + serialized_size(&pi.lv.iip_z.v_hat_tau_1) +
                     serialized_size(&pi.lv.nz.q0_tau_1) + serialized_size(&pi.lv.nz.w_tau_2) +
                     serialized_size(&pi.lv.p_tau_1) + serialized_size(&pi.lv.h_tau_1) +
                     serialized_size(&pi.lv.a_tau_1) + serialized_size(&pi.lv.b_tau_1) +
                     serialized_size(&pi.lv.c_tau_1) + serialized_size(&pi.lv.b_tau_2) +
                     serialized_size(&pi.lv.w_hat_tau_1);
    println!("LV Proof: {}", proof_size);

    // --- Encryptor's public LV params and header (no witness needed) ---
    let params = we::lv_public_linear_params(&crs, &dg.lv);
    let (hdr, key_enc) = we::lv_make_header(&params, &crs, &mut rng);
    
    // Header size: manually calculate
    let mut header_size = 0;
    for elem in &hdr.c1 {
        match elem {
            we::HeaderElem::G1(g) => header_size += serialized_size(g),
            we::HeaderElem::G2(g) => header_size += serialized_size(g),
        }
    }
    println!("Header: {}", header_size);

    // --- AEAD encrypt ---
    let mut msg = b"hello secret world".to_vec();
    let nonce: [u8; 12] = rng.random();
    
    let enc_start = Instant::now();
    let tag: Vec<u8> = aead_encrypt(&crs, &params.shape, &hdr, key_enc, nonce, &mut msg);
    let enc_time = enc_start.elapsed();
    
    let ciphertext_size = msg.len();
    let tag_size = tag.len();
    println!("Ciphertext: {}", ciphertext_size);
    println!("Tag: {}", tag_size);
    
    println!("\n=== TIMING ===");
    println!("Encryption: {:?}", enc_time);

    // --- Decryptor derives key from π + header, then decrypt ---
    let mut ct: Vec<u8> = msg.clone();
    
    let dec_start = Instant::now();
    let maybe_pt = decrypt_with_lv_header(&crs, &dg.lv, &params, &hdr, &pi.lv, nonce, &mut ct, &tag);
    let dec_time = dec_start.elapsed();
    
    println!("Decryption: {:?}", dec_time);
    
    match maybe_pt {
        Some(pt) => println!("\n=== RESULT ===\nDecryption OK: {}", String::from_utf8_lossy(&pt)),
        None => println!("\n=== RESULT ===\nDecryption failed"),
    }
}