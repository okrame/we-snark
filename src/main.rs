// src/main.rs
mod inner;
mod scs;
mod lv_gadgets;
mod we;
mod routing;

use ark_bw6_761::Fr as OuterFr;

fn main() -> anyhow::Result<()> {
    println!("╔════════════════════════════════════════╗");
    println!("║   WE CONSTRUCTION 2 SETUP              ║");
    println!("╚════════════════════════════════════════╝");
    
    let n = 8usize;
    let nmax = 2 * n;
    let (pp, crs) = we::setup(&mut rand::thread_rng(), n, nmax);
    println!("✓ SCS CRS generated (powers-of-τ)");
    println!("  Domain size: {}", n);
    println!("  Max degree: {}", nmax);
    
    let s_vals: Vec<OuterFr> = (0..n)
        .map(|i| OuterFr::from((42 + i) as u64))
        .collect();
    let (idx, _aux) = we::auxgen(&pp, &crs, &s_vals, nmax);
    println!("✓ AuxGen completed: IIP index with {} elements", s_vals.len());
    
    let (ed, dd) = we::digest(&pp, &crs, &idx);
    println!("✓ Digest computed: EncDigest + DecDigest");
    println!();

    println!("╔════════════════════════════════════════╗");
    println!("║   WITNESS ENCRYPTION                   ║");
    println!("╚════════════════════════════════════════╝");
    
    let msg = b"hello recursive WE with proper Garg construction";
    let ct = we::enc(&mut rand::thread_rng(), &ed, msg);
    println!("✓ Message encrypted");
    println!("  Plaintext: {:?}", String::from_utf8_lossy(msg));
    println!("  Ciphertext ct₁ slots: {}", ct.c1_g1.len());
    println!();

    println!("╔════════════════════════════════════════╗");
    println!("║   LV PROOF GENERATION                  ║");
    println!("╚════════════════════════════════════════╝");
    
    let w: Vec<OuterFr> = (0..n)
        .map(|i| OuterFr::from((i + 1) as u64))
        .collect();
    
    println!("→ Precomputing Lagrange basis...");
    let lag = lv_gadgets::precompute_lagrange(&crs.ck.domain);
    println!("✓ Lagrange basis precomputed");
    
    println!("→ Generating LV proof components...");
    let pi = we::prove(&dd, &idx, &crs.ck, &w, &lag)?;
    println!("✓ LV proof generated with components:");
    println!("  - [Q_X(τ)]₁  (indexed inner product quotient)");
    println!("  - [Q_Z(τ)]₁  (vanishing polynomial quotient)");
    println!("  - [w(τ)]₂    (witness commitment in G2)");
    println!("  - [v(τ)]₁    (inner product commitment in G1)");
    println!("  - Zero-check proof (NonZero gadget)");
    println!();

    println!("╔════════════════════════════════════════╗");
    println!("║   WITNESS DECRYPTION                   ║");
    println!("╚════════════════════════════════════════╝");
    
    println!("→ Performing linear verification...");
    let pt = we::dec(&dd, &ct, &pi)?;
    println!("✓ Linear verification equations satisfied:");
    println!("  e(C,w) = e(v,[y*⁻¹]) · e([Q_X],[τ-x*]) · e([Q_Z],Z)");
    println!("→ Key derived via H(s·b)");
    println!("→ AEAD decryption...");
    println!("✓ Decrypted plaintext: {:?}", String::from_utf8_lossy(&pt));
    
    assert_eq!(pt, msg, "Decryption mismatch!");
    println!("✓ Decryption successful - plaintexts match!");
    println!();

    Ok(())
}