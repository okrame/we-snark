mod scs;
mod iip;
mod nonzero;
mod verifier;
mod we;

use ark_bn254::Fr;
use rand::{rng, Rng};

use scs::CRS;
use iip::{iip_digest, iip_prove};
use nonzero::nonzero_prove;
use verifier::{LVDigest, LVProof};
use we::{aead_encrypt, decrypt_with_lv_header};

fn main() {
    let mut rng = rng();

    // --- Parameters ---
    // Domain size n = 4: slots [x, y, z, 1]
    let n = 4;
    let crs = CRS::setup(&mut rng, n);

    // Public index s (we use s = [0,0,1,0] so v = [z]_1; the multiplication relation is encoded
    // in the field-side construction of B(X) with (x, y, z) chosen s.t. x*y=z).
    let s = vec![Fr::from(0u32), Fr::from(0u32), Fr::from(1u32), Fr::from(0u32)];

    // Instance: choose (x, y), set z = x*y
    let x = Fr::from(12u32);
    let y = Fr::from(17u32);
    let z = x * y;

    // Witness evals on D: [x, y, z, 1]
    let w = vec![x, y, z, Fr::from(1u32)];

    // --- Digest / Prove ---
    let iip_vk = iip_digest(&crs, &s);
    let iip_pi = iip_prove(&crs, &s, &w);
    let nz_pi = nonzero_prove(&crs, &w, 3);

    let dg = LVDigest { iip: iip_vk, one_idx: 3 };
    let pi = LVProof { iip: iip_pi, nz: nz_pi };

    // --- Encryptor's public LV params and header (no witness needed) ---
    let params = we::lv_public_linear_params(&crs, &dg);
    let (hdr, key_enc) = we::lv_make_header(&params, &crs, &mut rng);

    println!("A_LV = {:?}", params.shape.a);
    println!("b_LV[3] = {:?}", params.shape.b[3]);

    // --- AEAD encrypt ---
    let mut msg = b"hello world".to_vec();
    let nonce: [u8; 12] = rng.random();
    let tag: Vec<u8> = aead_encrypt(&crs, &params.shape, &hdr, key_enc, nonce, &mut msg);


    // --- Decryptor derives key from π + header, then decrypt ---
    let mut ct: Vec<u8> = msg.clone(); // after aead.enc, msg is mutated to the ciphertext. so here we copy its bytes, not the original message.
    let maybe_pt = decrypt_with_lv_header(&crs, &dg, &params, &hdr, &pi, nonce, &mut ct, &tag);
    match maybe_pt {
        Some(pt) => println!("Decryption OK: {}", String::from_utf8_lossy(&pt)),
        None => println!("Decryption failed"),
    }
}