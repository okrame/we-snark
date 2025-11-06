mod circuits;
mod groth16_lv;
mod types;
mod wit_dec;
mod wit_enc;

use ark_bn254::{Bn254, Fr};
use ark_groth16::Groth16;
use ark_snark::SNARK;
use circuits::simple_mul::MulCircuit;
use rand::thread_rng;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = thread_rng();

    // Setup
    let circuit = MulCircuit::<Fr> {
        x: None,
        y: None,
        z: None,
        lambda1: None,
        lambda2: None,
        lambda3: None,
    };
    
    let (pk, vk) = Groth16::<Bn254>::circuit_specific_setup(circuit, &mut rng)?;
    
    println!("╔════════════════════════════════════════╗");
    println!("║        CIRCUIT SETUP                   ║");
    println!("╚════════════════════════════════════════╝");
    println!("Number of constraints: {}", pk.vk.gamma_abc_g1.len() - 1);
    println!();

    // Witness
    let x = Fr::from(3u64);
    let y = Fr::from(4u64);

    // Generate LV proof
    println!("╔════════════════════════════════════════╗");
    println!("║        PROOF GENERATION                ║");
    println!("╚════════════════════════════════════════╝");
    let (lv_proof, z) = groth16_lv::lv_prove_mul(&pk, x, y)?;
    println!(" LV proof generated for z = {}", z);
    println!();

    // Encrypt a message
    println!("╔════════════════════════════════════════╗");
    println!("║        ENCRYPTION                      ║");
    println!("╚════════════════════════════════════════╝");
    let msg = b"my secret message";
    let ct = wit_enc::we_encrypt(&vk, z, msg)?;
    println!();

    // Decrypt using the proof
    println!("╔════════════════════════════════════════╗");
    println!("║        DECRYPTION                      ║");
    println!("╚════════════════════════════════════════╝");
    let plaintext = wit_dec::we_decrypt(&vk, &ct, &lv_proof)?;
    assert_eq!(plaintext, msg);
    println!();
    println!("Message recovered: {:?}", String::from_utf8_lossy(&plaintext));

    Ok(())
}