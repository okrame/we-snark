mod circuits;
mod inner;
mod recursion;
mod lv;
mod we_lv;

use ark_bls12_377::Fr as InnerFr;
use circuits::simple_mul::MulCircuit;

fn main() -> anyhow::Result<()> {
    println!("╔════════════════════════════════════════╗");
    println!("║   RECURSIVE SNARK-BASED WE             ║");
    println!("╚════════════════════════════════════════╝");
    println!();

    println!("╔════════════════════════════════════════╗");
    println!("║   INNER SETUP (BLS12-377)              ║");
    println!("╚════════════════════════════════════════╝");
    
    let x = InnerFr::from(3u64);
    let y = InnerFr::from(4u64);
    let z = x * y;
    
    let dummy_circuit = MulCircuit::<InnerFr> {
        x: None,
        y: None,
        z: None,
    };
    
    let (pk_in, vk_in) = inner::setup(dummy_circuit)?;
    println!("Inner Groth16 keys generated (BLS12-377)");
    println!("Instance: z = {}", z);
    println!();

    println!("╔════════════════════════════════════════╗");
    println!("║   INNER PROOF GENERATION               ║");
    println!("╚════════════════════════════════════════╝");
    
    let circuit = MulCircuit {
        x: Some(x),
        y: Some(y),
        z: Some(z),
    };
    
    let inner_proof = inner::prove(&pk_in, circuit)?;
    println!("Inner proof generated for x={}, y={}, z={}", x, y, z);
    println!();

    println!("╔════════════════════════════════════════╗");
    println!("║   OUTER SETUP (BW6-761)                ║");
    println!("╚════════════════════════════════════════╝");
    
    let dummy_outer = circuits::outer_verify_inner::OuterVerifiesInnerCircuit {
        inner_vk: vk_in.clone(),
        inner_public_inputs: vec![InnerFr::from(1u64)],
        inner_proof: inner_proof.clone(),
    };
    
    let (pk_out, vk_out) = recursion::setup_outer(dummy_outer)?;
    println!("Outer Groth16 keys generated (BW6-761)");
    println!("Outer circuit verifies inner proofs");
    println!();

    println!("╔════════════════════════════════════════╗");
    println!("║   OUTER PROOF (RECURSIVE)              ║");
    println!("╚════════════════════════════════════════╝");
    
    let outer_proof = recursion::prove_outer(
        &pk_out,
        vk_in.clone(),
        vec![z],
        inner_proof.clone(),
    )?;
    println!("Outer proof generated");
    println!("Proves: 'inner proof for z={} is valid'", z);
    println!();

    println!("╔════════════════════════════════════════╗");
    println!("║   LV EXTRACTION                        ║");
    println!("╚════════════════════════════════════════╝");
    
    let lin = lv::derive_a_from_outer_proof(&vk_out, &outer_proof, 4);
    println!("Linear verification material extracted");
    println!("LV slots: {}", lin.y_slots.len());
    println!();

    println!("╔════════════════════════════════════════╗");
    println!("║   WITNESS ENCRYPTION                   ║");
    println!("╚════════════════════════════════════════╝");
    
    let ed = we_lv::setup_lv_slots(4);
    let msg = b"hello recursive WE";
    let ct = we_lv::enc(&ed, msg);
    println!("Message encrypted: {:?}", String::from_utf8_lossy(msg));
    println!("Ciphertext slots: {}", ct.s_g1.len());
    println!();

    println!("╔════════════════════════════════════════╗");
    println!("║   WITNESS DECRYPTION                   ║");
    println!("╚════════════════════════════════════════╝");
    
    let pt = we_lv::dec(&ct, &lin);
    println!("Decrypted: {:?}", String::from_utf8_lossy(&pt));
    assert_eq!(pt, msg);
    println!("✓ Decryption successful!");
    println!();

    Ok(())
}