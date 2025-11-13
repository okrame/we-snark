//src/verifier.rs
use crate::iip::{IIPDigest, IIPProof, iip_verify};
use crate::nonzero::{NonZeroProof, nonzero_verify};
use crate::scs::CRS;
use ark_bn254::{Bn254, Fq12, Fr, G1Projective as G1, G2Projective as G2};
use ark_ec::pairing::Pairing;
use ark_ec::PrimeGroup;
use ark_ff::Field;
use ark_ff::One;
use ark_ff::PrimeField;
use ark_poly::EvaluationDomain;

#[derive(Clone, Copy)]
pub enum ColSide { ProofG1PublicG2, ProofG2PublicG1 }

#[derive(Clone)]
pub struct LVColMeta {
    pub side: ColSide,
    pub g1_pub: Option<G1>,
    pub g2_pub: Option<G2>,
}

pub enum ProofElem { G1(G1), G2(G2) }

#[derive(Clone)] 
pub struct LVDigest {
    pub iip: IIPDigest,
    pub one_idx: usize,
}

pub struct LVCoords(pub [Fq12; LV_NUM_COORDS]);
pub(crate) fn build_lv_coords(crs: &CRS, dg: &LVDigest, pi: &LVProof) -> Option<LVCoords> {
    // The NonZero and IIP commitments to B(τ) must match
    if pi.iip.w_tau_2 != pi.nz.w_tau_2 { return None; }

    let g1 = <Bn254 as Pairing>::G1::generator();
    let g2 = <Bn254 as Pairing>::G2::generator();

    // y*^{-1}
    let y_inv = dg.iip.y_star.inverse().unwrap();

    // d = D[one_idx]; [τ - d]_2
    let d = crs.domain.element(dg.one_idx);
    let tau_minus_d_2 = crs.g2_tau_pow(1) - g2.mul_bigint(d.into_bigint());

    // Fill the 10 coordinates (PairingOutputs turned into Fq12)
    let c0 = <Bn254 as Pairing>::pairing(dg.iip.C,                pi.iip.w_tau_2).0;
    let c1 = <Bn254 as Pairing>::pairing(pi.iip.v_g1.mul_bigint(y_inv.into_bigint()), g2).0;
    let c2 = <Bn254 as Pairing>::pairing(pi.iip.QX_tau_1,         dg.iip.tau_2).0;
    let c3 = <Bn254 as Pairing>::pairing(pi.iip.QZ_tau_1,         dg.iip.Z_tau_2).0;
    let c4 = <Bn254 as Pairing>::pairing(pi.iip.QX_tau_1,         dg.iip.tau_N_minus_n_plus_2_2).0;
    let c5 = <Bn254 as Pairing>::pairing(pi.iip.QX_hat_tau_1,     g2).0;
    let c6 = <Bn254 as Pairing>::pairing(pi.iip.v_g1,             dg.iip.tau_N_2).0;
    let c7 = <Bn254 as Pairing>::pairing(pi.iip.v_hat_tau_1,      g2).0;
    let c8 = <Bn254 as Pairing>::pairing(g1,                      pi.nz.w_tau_2).0;
    let c9 = <Bn254 as Pairing>::pairing(pi.nz.q0_tau_1,          tau_minus_d_2).0;

    Some(LVCoords([c0,c1,c2,c3,c4,c5,c6,c7,c8,c9]))
}

/// Collect proof-side elements per column (G1 or G2), matching column order
pub(crate) fn build_proof_side_elems(_crs: &CRS, dg: &LVDigest, pi: &LVProof)
    -> Option<[ProofElem; LV_NUM_COORDS]>
{
    if pi.iip.w_tau_2 != pi.nz.w_tau_2 { return None; }

    let y_inv = dg.iip.y_star.inverse().unwrap();

    Some([
        ProofElem::G2(pi.iip.w_tau_2),
        ProofElem::G1(pi.iip.v_g1.mul_bigint(y_inv.into_bigint())),
        ProofElem::G1(pi.iip.QX_tau_1),
        ProofElem::G1(pi.iip.QZ_tau_1),
        ProofElem::G1(pi.iip.QX_tau_1),
        ProofElem::G1(pi.iip.QX_hat_tau_1),
        ProofElem::G1(pi.iip.v_g1),
        ProofElem::G1(pi.iip.v_hat_tau_1),
        ProofElem::G2(pi.nz.w_tau_2),
        ProofElem::G1(pi.nz.q0_tau_1),
    ])
}

#[derive(Clone)]
pub struct LVProof {
    pub iip: IIPProof,
    pub nz: NonZeroProof,
    pub w: Vec<Fr>,
}

/// Number of GT-coordinates we use in A_LV · π = b_LV.
pub const LV_NUM_COORDS: usize = 10;

/// A_LV and b_LV as described above.
/// - a[i][j] ∈ {-1,0,1} describes exponent α_{i,j} on coordinate c_j in equation i.
/// - b[i] ∈ GT is the RHS constant for equation i.
pub struct LVShape {
    pub rows: usize,
    pub a: [[i8; LV_NUM_COORDS]; 4], // here rows=4 fixed, you can generalize later
    pub b: [Fq12; 4],
}

impl LVDigest {
    /// Build the A_LV and b_LV used by the LV verifier, in the abstract
    /// GT-coordinate basis c_0..c_9 
    pub fn linear_shape(&self, _crs: &CRS) -> LVShape {
        let rows = 4;

        // Matrix A_LV: fill explicitly with the α_{i,j} above.
        let mut a = [[0i8; LV_NUM_COORDS]; 4];

        // Eq 0: c0 * c1^{-1} * c2^{-1} * c3^{-1} = 1
        a[0] = [ 1, -1, -1, -1,  0,  0,  0,  0,  0,  0];

        // Eq 1: c4 * c5^{-1} = 1
        a[1] = [ 0,  0,  0,  0,  1, -1,  0,  0,  0,  0];

        // Eq 2: c6 * c7^{-1} = 1
        a[2] = [ 0,  0,  0,  0,  0,  0,  1, -1,  0,  0];

        // Eq 3: c8 * c9^{-1} = e(g1,g2)
        a[3] = [ 0,  0,  0,  0,  0,  0,  0,  0,  1, -1];

        // b_LV: all 1 except the last, which is e(g1,g2)
        let gt_one = Fq12::one();
        let gt_const: Fq12 = <Bn254 as Pairing>::pairing(
            <Bn254 as Pairing>::G1::generator(),
            <Bn254 as Pairing>::G2::generator(),
        ).0;

        let b = [gt_one.clone(), gt_one.clone(), gt_one.clone(), gt_const];

        LVShape { rows, a, b }
    }

    /// Map each column to its public base and orientation
    pub fn column_metadata(&self, crs: &CRS) -> [LVColMeta; LV_NUM_COORDS] {
        let g1 = <Bn254 as Pairing>::G1::generator();
        let g2 = <Bn254 as Pairing>::G2::generator();
        let d = crs.domain.element(self.one_idx);
        let tau_minus_d_2 = crs.g2_tau_pow(1) - g2.mul_bigint(d.into_bigint());

        [
            // c0 = e(C, w_tau_2): proof is G2, public base is G1 (C)
            LVColMeta { side: ColSide::ProofG2PublicG1, g1_pub: Some(self.iip.C), g2_pub: None },
            // c1 = e(v_g1 * y_inv, g2): proof is G1, public base is g2
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(g2) },
            // c2 = e(QX_tau_1, tau_2): proof G1, public G2
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(self.iip.tau_2) },
            // c3 = e(QZ_tau_1, Z_tau_2): proof G1, public G2
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(self.iip.Z_tau_2) },
            // c4 = e(QX_tau_1, tau_{N-n+2,2})
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(self.iip.tau_N_minus_n_plus_2_2) },
            // c5 = e(QX_hat_tau_1, g2)
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(g2) },
            // c6 = e(v_g1, tau_N_2)
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(self.iip.tau_N_2) },
            // c7 = e(v_hat_tau_1, g2)
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(g2) },
            // c8 = e(g1, w_tau_2): proof G2, public G1
            LVColMeta { side: ColSide::ProofG2PublicG1, g1_pub: Some(g1), g2_pub: None },
            // c9 = e(q0_tau_1, (tau - d)_2)
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(tau_minus_d_2) },
        ]
    }
}


pub fn recover_sb_via_linear_check(
    shape: &LVShape,
    coords: &[Fq12; LV_NUM_COORDS],
) -> bool {
    for i in 0..shape.rows {
        let mut lhs = Fq12::one();
        for j in 0..LV_NUM_COORDS {
            let e = shape.a[i][j];
            if e == 0 { continue; }
            if e == 1  { lhs *= &coords[j]; }
            if e == -1 {
                let inv = coords[j].inverse().unwrap();
                lhs *= &inv;
            }
        }
        if lhs != shape.b[i] { return false; }
    }
    true
}

pub fn lv_verify(crs: &CRS, dg: &LVDigest, pi: &LVProof) -> bool {
    // Basic relation checks on the explicit witness vector.
    if pi.w.len() != 4 {
        return false;
    }
    let x   = pi.w[0];
    let y   = pi.w[1];
    let z   = pi.w[2];
    let one = pi.w[3];

    //Enforce the NP relation x * y = z and w_3 = 1 in Fr.
    if x * y != z {
        return false;
    }

    // Non-zero slot: w_3 must literally be 1
    if !one.is_one() {
        return false;
    }

    // Optional: keep the original gadgets as safety checks in debug builds
    #[cfg(debug_assertions)]
    {
        if !iip_verify(&dg.iip, &pi.iip) { return false; }
        if !nonzero_verify(crs, &pi.nz, dg.one_idx) { return false; }
    }

    let shape = dg.linear_shape(crs);
    let coords = match build_lv_coords(crs, dg, pi) {
        Some(c) => c,
        None => return false,
    };

    recover_sb_via_linear_check(&shape, &coords.0)
}
