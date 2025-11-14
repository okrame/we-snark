# Mini Witness Encryption with LV-SNARK

This repo is a tiny end-to-end prototype of Witness Encryption built from a linearly verifiable (LV) SNARK-like gadget stack (see this https://eprint.iacr.org/2025/1364).
The NP relation supported is:

$$
x \cdot y = z.
$$

Everything runs on BN254 using arkworks crates.


> Note: This captures a tiny linearly verifiable gadget stack for a fixed constraint system and should not be confused with a general-purpose WE construction. Standard bilinear groups only provide degree-2 structure, so collapsing the verification of arbitrary circuits into a single GT-linear equation would require stronger primitives such as multilinear maps.

## What’s going on?

We build a set of GT-linear gadgets (IIP, NonZero, Mul, QAP-binding) whose combined verification collapses into a single linear check:

$$
A_{\text{LV}} \cdot \pi = b_{\text{LV}}.
$$

Here $\pi$ is a vector of GT coordinates derived from the witness.

* The encryptor uses only the public LV digest $(A_{\text{LV}}, b_{\text{LV}})$ to encrypt.
  → Encryptor never sees $(x,y,z)$ and never constructs a proof.

* The decryptor knows the witness $(x,y,z)$, locally builds the LV proof $\pi$, checks the LV system, and recovers the symmetric key.

If you don’t know a valid $(x,y,z)$ with $x \cdot y = z$, you cannot reconstruct $\pi$, and the key stays hidden.

---

## Gadget stack

* **IIP gadgets** (selectors for $x$, $y$, $z$):
  expose the three witness values in GT linearly.

* **NonZero gadget**:
  enforces that the last slot of the witness is $1$ (classic trick).

* **Mul gadget**:
  enforces multiplicative consistency using a QAP-style constraint:

  * $A(\tau) = x$
  * $B(\tau) = y$
  * $C(\tau) = z$
  * $P(\tau) = A(\tau) B(\tau) - C(\tau)$
  * $P(\tau) = H(\tau) Z(\tau)$

* **C–z binding**:
  ensures that the polynomial output $C$ encodes the same $z$ exposed by the IIP gadget.

All these constraints are flattened into GT coordinates `c0 … c15` and a handful of LV rows.

No field-side `x * y == z` check is needed — the gadgets enforce it.

---

## Encryption flow

Encryptor only needs:

* the LV digest (matrix shape + CRS-derived GT bases),
* a message.

**Steps**

1. Sample randomness $s$
2. Compute
   $$c_1 = s \cdot A_{\text{LV}},\qquad
   c_2 = s \cdot b_{\text{LV}} + k$$
   where $k$ is the symmetric key.
3. Encrypt the plaintext using AES-GCM with key $k$.
4. Output:

   ```
   { LVHeader, c1, c2, AES_nonce, AES_tag, ciphertext }
   ```

---

## Decryption flow

Decryptor knows the witness $(x, y, z)$.

**Steps**

1. Build the witness vector
   $$w = [x, y, z, 1]$$
2. Run the gadgets to compute the LV proof $\pi$.
3. Check linear consistency:
   $$A_{\text{LV}} \cdot \pi = b_{\text{LV}}$$
4. Recover
   $$k = c_1 \cdot \pi - c_2$$
5. AES-GCM decrypt.

If the witness is invalid, the LV system fails and the key cannot be reconstructed.

## Language-level vs Instance-level WE

Right now, the LV digest represents the **language**

$$L = {(x,y,z) \mid x \cdot y = z},$$

so the encryptor does not need to know any specific $z_0$.

To obtain *instance-specific* WE (e.g., “decrypt only if $x \cdot y = z_0$”), we can extend the LV digest to bake in the public value $z_0$.
The encryptor still remains oblivious to $(x,y)$ and never handles a witness.
