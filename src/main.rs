mod circuits;
mod groth16_lv;

use ark_bn254::{Bn254, Fr};
use ark_groth16::{Groth16, Proof};
use ark_snark::SNARK;
use circuits::simple_mul::MulCircuit;
use rand::thread_rng;

use groth16_lv::{lv_prove_mul, lv_instance_from_vk_and_z, lv_verify_mul};


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = thread_rng();

    let params = Groth16::<Bn254>::circuit_specific_setup(
        MulCircuit::<Fr> { x: None, y: None },
        &mut rng,
    )?;

    let x = Fr::from(3u64);
    let y = Fr::from(4u64);
    let z = x * y; 

    let circuit = MulCircuit {
        x: Some(x),
        y: Some(y),
    };

    let proof: Proof<Bn254> = Groth16::<Bn254>::prove(&params.0, circuit, &mut rng)?;

    // verifying key
    let pvk = Groth16::<Bn254>::process_vk(&params.1)?;

    // verify the proof against the public input z
    let verified = Groth16::<Bn254>::verify_with_processed_vk(&pvk, &[z], &proof)?;
    assert!(verified);

    println!("Proof verified successfully!");

    // extract the linear pairing system for this (vk, u, π)
    let vk = &params.1;

    // --- LV PROVE ---
    let (lv_proof, z) = lv_prove_mul(&params.0, x, y)?;
    println!("z = x * y = {}", z);

    // --- LV INSTANCE SIDE ---
    let inst = lv_instance_from_vk_and_z(&vk, z);

    // --- LV VERIFY ---
    let ok = lv_verify_mul(&vk, &inst, &lv_proof);
    println!("LV-SNARK verify: {}", ok);


    Ok(())
}