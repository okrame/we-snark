mod scs;
mod iip;
mod nonzero;
mod verifier;
mod we;
mod mul_snark;

use ark_bn254::Fr;
use rand::{rng, Rng};

use scs::CRS;
use we::{aead_encrypt, decrypt_with_lv_header};
use mul_snark::{MulDigest, MulWitness, mul_prove, mul_verify};

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

    let dg = MulDigest::setup(&crs);
    let pi = mul_prove(&crs, &dg, &w);

    // sanity check
    assert!(mul_verify(&crs, &dg, &pi));

    // --- Encryptor's public LV params and header (no witness needed) ---
    let params = we::lv_public_linear_params(&crs, &dg.lv);
    let (hdr, key_enc) = we::lv_make_header(&params, &crs, &mut rng);

    println!("A_LV = {:?}", params.shape.a);
    println!("b_LV[3] = {:?}", params.shape.b[3]);

    // --- AEAD encrypt ---
    let mut msg = b"hello world".to_vec();
    let nonce: [u8; 12] = rng.random();
    let tag: Vec<u8> = aead_encrypt(&crs, &params.shape, &hdr, key_enc, nonce, &mut msg);


    // --- Decryptor derives key from π + header, then decrypt ---
    let mut ct: Vec<u8> = msg.clone(); // after aead.enc, msg now holds ciphertext bytes
    let maybe_pt = decrypt_with_lv_header(&crs, &dg.lv, &params, &hdr, &pi.lv, nonce, &mut ct, &tag);
    match maybe_pt {
        Some(pt) => println!("Decryption OK: {}", String::from_utf8_lossy(&pt)),
        None => println!("Decryption failed"),
    }
}