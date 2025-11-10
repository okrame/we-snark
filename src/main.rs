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
use we::{aead_encrypt, decrypt_with_lv, derive_key_from_lv};

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
    let nz_pi = nonzero_prove(&crs, &w, 3); // enforce w[3] == 1

    let dg = LVDigest { iip: iip_vk, one_idx: 3 };
    let pi = LVProof { iip: iip_pi, nz: nz_pi };

    // --- Derive key from LV (what decryptor will do) ---
    let key = derive_key_from_lv(&dg.iip, &pi.iip);
    let mut msg = b"hello, LV world".to_vec();
    let nonce: [u8;12] = rng.random();

    // Encrypt (encryptor does NOT know the witness; here we just reuse the same derived key
    // to demonstrate end-to-end; in a full WE, the encryptor uses the LV verifier’s linear
    // coefficients to produce a header that makes this key recoverable).
    let tag: Vec<u8> = aead_encrypt(key, nonce, &mut msg, b"AAD");
    
    let lv_params = we::lv_public_linear_params(&crs, &dg);
    println!("A_LV = {:?}", lv_params.shape.a);
    println!("b_LV[3] = {:?}", lv_params.shape.b[3]);

    // --- Decrypt using witness (x,y) that satisfies x*y=z ---
    let mut ct = msg.clone();
    let maybe_pt = decrypt_with_lv(&dg, &pi, nonce, &mut ct, &tag, b"AAD");
    match maybe_pt {
        Some(pt) => println!("Decryption OK: {}", String::from_utf8_lossy(&pt)),
        None => println!("Decryption failed"),
    }
}