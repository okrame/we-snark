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
#[allow(non_snake_case)]
pub struct LVDigest {
    pub iip_x: IIPDigest, 
    pub iip_y: IIPDigest, 
    pub iip_z: IIPDigest,
    pub one_idx: usize,
    pub mul_z_tau_2: G2,
    pub instance_z: Fr,
    // MaxDeg parameters for the IIP witness polynomial B(X)
    pub d_bound: usize,     // e.g. n-1
    pub tau_N_minus_d_1: G1 // [τ^{N-d}]_1
}

pub struct LVCoords(pub [Fq12; LV_NUM_COORDS]);
pub(crate) fn build_lv_coords(crs: &CRS, dg: &LVDigest, pi: &LVProof) -> Option<LVCoords> {
    // The NonZero and IIP commitments to B(τ) must match
    if pi.iip_z.w_tau_2 != pi.nz.w_tau_2 { return None; }

    let g1 = <Bn254 as Pairing>::G1::generator();
    let g2 = <Bn254 as Pairing>::G2::generator();

    // y*^{-1}
    let y_inv = dg.iip_z.y_star.inverse().unwrap();

    // d = D[one_idx]; [τ - d]_2
    let d = crs.domain.element(dg.one_idx);
    let tau_minus_d_2 = crs.g2_tau_pow(1) - g2.mul_bigint(d.into_bigint());

    // Fill the coordinates (PairingOutputs turned into Fq12)
    let c0 = <Bn254 as Pairing>::pairing(dg.iip_z.C,                pi.iip_z.w_tau_2).0;
    let c1 = <Bn254 as Pairing>::pairing(pi.iip_z.v_g1.mul_bigint(y_inv.into_bigint()), g2).0;
    let c2 = <Bn254 as Pairing>::pairing(pi.iip_z.QX_tau_1,         dg.iip_z.tau_2).0;
    let c3 = <Bn254 as Pairing>::pairing(pi.iip_z.QZ_tau_1,         dg.iip_z.Z_tau_2).0;
    let c4 = <Bn254 as Pairing>::pairing(pi.iip_z.QX_tau_1,         dg.iip_z.tau_N_minus_n_plus_2_2).0;
    let c5 = <Bn254 as Pairing>::pairing(pi.iip_z.QX_hat_tau_1,     g2).0;
    let c6 = <Bn254 as Pairing>::pairing(pi.iip_z.v_g1,             dg.iip_z.tau_N_2).0;
    let c7 = <Bn254 as Pairing>::pairing(pi.iip_z.v_hat_tau_1,      g2).0;
    let c8 = <Bn254 as Pairing>::pairing(g1,                      pi.nz.w_tau_2).0;
    let c9 = <Bn254 as Pairing>::pairing(pi.nz.q0_tau_1,          tau_minus_d_2).0;
    // Mul-gadget coordinates
    let c10 = <Bn254 as Pairing>::pairing(pi.p_tau_1, g2).0;
    let c11 = <Bn254 as Pairing>::pairing(pi.h_tau_1, dg.mul_z_tau_2).0;
    let c12 = <Bn254 as Pairing>::pairing(pi.a_tau_1, g2).0; 
    let c13 = <Bn254 as Pairing>::pairing(pi.b_tau_1, g2).0; 

    // C–z binding coordinates:
    // c14 = e(v_g1, g2), where v_g1 = z from IIP selector s = [0,0,1,0]
    // c15 = e(C(τ)_1, g2), where C(X) = z is the QAP output polynomial
    let c14 = <Bn254 as Pairing>::pairing(pi.iip_z.v_g1, g2).0;
    let c15 = <Bn254 as Pairing>::pairing(pi.c_tau_1, g2).0;

    // --- MaxDeg gadget coordinates ---
    // c16 = e([τ^{N-d}]_1, [B(τ)]_2) where B(X) is the IIP witness polynomial
    let c16 = <Bn254 as Pairing>::pairing(dg.tau_N_minus_d_1, pi.iip_z.w_tau_2).0;
    // c17 = e([X^{N-d} B(X)]_1, g2)
    let c17 = <Bn254 as Pairing>::pairing(pi.w_hat_tau_1, g2).0;

    Some(LVCoords([
    c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,
    c10,c11,c12,c13,c14,c15,c16,c17
]))
}

/// Collect proof-side elements per column (G1 or G2), matching column order
pub(crate) fn build_proof_side_elems(_crs: &CRS, dg: &LVDigest, pi: &LVProof)
    -> Option<[ProofElem; LV_NUM_COORDS]>
{
    if pi.iip_z.w_tau_2 != pi.nz.w_tau_2 { return None; }

    let y_inv = dg.iip_z.y_star.inverse().unwrap();

    Some([
        ProofElem::G2(pi.iip_z.w_tau_2),
        ProofElem::G1(pi.iip_z.v_g1.mul_bigint(y_inv.into_bigint())),
        ProofElem::G1(pi.iip_z.QX_tau_1),
        ProofElem::G1(pi.iip_z.QZ_tau_1),
        ProofElem::G1(pi.iip_z.QX_tau_1),
        ProofElem::G1(pi.iip_z.QX_hat_tau_1),
        ProofElem::G1(pi.iip_z.v_g1),
        ProofElem::G1(pi.iip_z.v_hat_tau_1),
        ProofElem::G2(pi.nz.w_tau_2),
        ProofElem::G1(pi.nz.q0_tau_1),
        // Mul gadget (P,H,A,B)
        ProofElem::G1(pi.p_tau_1),
        ProofElem::G1(pi.h_tau_1),
        ProofElem::G1(pi.a_tau_1),
        ProofElem::G1(pi.b_tau_1),
        // C–z binding reuses v_z_g1 and C(τ)_1
        ProofElem::G1(pi.iip_z.v_g1),
        ProofElem::G1(pi.c_tau_1),
        // MaxDeg: witness B(τ) and shifted commitment
        ProofElem::G2(pi.iip_z.w_tau_2), // c16 proof element (matches ProofG2PublicG1)
        ProofElem::G1(pi.w_hat_tau_1),   // c17 proof element (matches ProofG1PublicG2)
    ])
}

#[derive(Clone)]
pub struct LVProof {
    pub iip_x: IIPProof,
    pub iip_y: IIPProof,
    pub iip_z: IIPProof,
    pub nz: NonZeroProof,
    pub w: Vec<Fr>,
    // Mul-gadget commitments
    pub p_tau_1: G1, // [P(τ)]_1
    pub h_tau_1: G1, // [H(τ)]_1
    pub a_tau_1: G1, // [A(τ)]_1
    pub b_tau_1: G1, // [B(τ)]_1  (for A/B binding)
    pub c_tau_1: G1, // [C(τ)]_1
    pub b_tau_2: G2, // [B(τ)]_2 for the P = A·B - C check
    pub w_hat_tau_1: G1,
}

/// Number of GT-coordinates we use in A_LV · π = b_LV.
pub const LV_NUM_COORDS: usize = 18;

/// A_LV and b_LV as described above.
/// - a[i][j] ∈ {-1,0,1} describes exponent α_{i,j} on coordinate c_j in equation i.
/// - b[i] ∈ GT is the RHS constant for equation i.
pub struct LVShape {
    pub rows: usize,
    pub a: [[i8; LV_NUM_COORDS]; 8], // here rows fixed, you can generalize later
    pub b: [Fq12; 8],
}

impl LVDigest {
        pub fn linear_shape(&self, _crs: &CRS) -> LVShape {
        let rows = 8;

        let mut a = [[0i8; LV_NUM_COORDS]; 8];

        // Eq 0: c0 * c1^{-1} * c2^{-1} * c3^{-1} = 1
        a[0] = [ 1, -1, -1, -1,  0,  0,  0,  0,  0,  0,
                 0,  0,  0,  0,  0,  0,  0,  0];

        // Eq 1: c4 * c5^{-1} = 1
        a[1] = [ 0,  0,  0,  0,  1, -1,  0,  0,  0,  0,
                 0,  0,  0,  0,  0,  0,  0,  0];

        // Eq 2: c6 * c7^{-1} = 1
        a[2] = [ 0,  0,  0,  0,  0,  0,  1, -1,  0,  0,
                 0,  0,  0,  0,  0,  0,  0,  0];

        // Eq 3: c8 * c9^{-1} = e(g1,g2)
        a[3] = [ 0,  0,  0,  0,  0,  0,  0,  0,  1, -1,
                 0,  0,  0,  0,  0,  0,  0,  0];

        // Eq 4 (Mul QAP): c10 * c11^{-1} = 1
        a[4] = [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                 1, -1,  0,  0,  0,  0,  0,  0];

        // Eq 5 (C–z binding): c14 * c15^{-1} = 1
        a[5] = [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                 0,  0,  0,  0,  1, -1,  0,  0];

        // Eq 6 (MaxDeg for B): c16 * c17^{-1} = 1
        a[6] = [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                 0,  0,  0,  0,  0,  0,  1, -1];

        // Eq 7 instance binding z = z0
        a[7] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 1, 0, 0, 0];


        let gt_one = Fq12::one();
        let gt_const: Fq12 = <Bn254 as Pairing>::pairing(
            <Bn254 as Pairing>::G1::generator(),
            <Bn254 as Pairing>::G2::generator(),
        ).0;

        let mut b = [gt_one.clone(); 8];
        b[3] = gt_const;

        // Eq 7: z = z0 ⇒ c14 = e(z0·G1, G2)
        let g1 = <Bn254 as Pairing>::G1::generator();
        let g2 = <Bn254 as Pairing>::G2::generator();
        let z0_g1 = g1.mul_bigint(self.instance_z.into_bigint());
        b[7] = <Bn254 as Pairing>::pairing(z0_g1, g2).0;
        
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
            LVColMeta { side: ColSide::ProofG2PublicG1, g1_pub: Some(self.iip_z.C), g2_pub: None },
            // c1 = e(v_g1 * y_inv, g2): proof is G1, public base is g2
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(g2) },
            // c2 = e(QX_tau_1, tau_2): proof G1, public G2
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(self.iip_z.tau_2) },
            // c3 = e(QZ_tau_1, Z_tau_2): proof G1, public G2
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(self.iip_z.Z_tau_2) },
            // c4 = e(QX_tau_1, tau_{N-n+2,2})
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(self.iip_z.tau_N_minus_n_plus_2_2) },
            // c5 = e(QX_hat_tau_1, g2)
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(g2) },
            // c6 = e(v_g1, tau_N_2)
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(self.iip_z.tau_N_2) },
            // c7 = e(v_hat_tau_1, g2)
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(g2) },
            // c8 = e(g1, w_tau_2): proof G2, public G1
            LVColMeta { side: ColSide::ProofG2PublicG1, g1_pub: Some(g1), g2_pub: None },
            // c9 = e(q0_tau_1, (tau - d)_2)
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(tau_minus_d_2) },
            // c10 = e(P_tau_1, g2)
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(g2) },
            // c11 = e(H_tau_1, Z_tau_2)
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(self.mul_z_tau_2) },
            // c12 = e(A_tau_1, g2) optional
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(g2) },
            // c13 = e(C_tau_1, g2) optional
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(g2) },
            // c14 = e(v_g1, g2)  (z from IIP)
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(g2) },
            // c15 = e(C_tau_1, g2)  (z from QAP C)
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(g2) },

            // c16 = e([τ^{N-d}]_1, [B(τ)]_2): proof G2, public G1
            LVColMeta { side: ColSide::ProofG2PublicG1, g1_pub: Some(self.tau_N_minus_d_1), g2_pub: None },

            // c17 = e([X^{N-d} B(X)]_1, g2): proof G1, public g2
            LVColMeta { side: ColSide::ProofG1PublicG2, g1_pub: None, g2_pub: Some(g2) },
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

#[allow(non_snake_case)]
pub fn lv_verify(crs: &CRS, dg: &LVDigest, pi: &LVProof) -> bool {
    // Basic relation check on witness length.
    if pi.w.len() != 4 {
        return false;
    }

    // Optional: keep the original gadgets as safety checks in debug builds
    #[cfg(debug_assertions)]
    {
        if !iip_verify(&dg.iip_x, &pi.iip_x) { return false; }
        if !iip_verify(&dg.iip_y, &pi.iip_y) { return false; }
        if !iip_verify(&dg.iip_z, &pi.iip_z) { return false; }
        if !nonzero_verify(crs, &pi.nz, dg.one_idx) { return false; }
    }

    // Optional: keep the original gadgets as safety checks in debug builds
    #[cfg(debug_assertions)]
    {
        if !iip_verify(&dg.iip_x, &pi.iip_x) { return false; }
        if !iip_verify(&dg.iip_y, &pi.iip_y) { return false; }
        if !iip_verify(&dg.iip_z, &pi.iip_z) { return false; }
        if !nonzero_verify(crs, &pi.nz, dg.one_idx) { return false; }
    }

    // group-level Mul gadget at τ:
    // A/B binding: enforce that A(τ) = x and B(τ) = y using IIP selectors
    {
        let g2 = <Bn254 as Pairing>::G2::generator();

        // From IIP with s_x = [1,0,0,0], v_x = x · [1]_1
        let gt_A = <Bn254 as Pairing>::pairing(pi.a_tau_1, g2).0;
        let gt_x = <Bn254 as Pairing>::pairing(pi.iip_x.v_g1, g2).0;
        if gt_A != gt_x {
            return false;
        }

        // From IIP with s_y = [0,1,0,0], v_y = y · [1]_1
        let gt_B = <Bn254 as Pairing>::pairing(pi.b_tau_1, g2).0;
        let gt_y = <Bn254 as Pairing>::pairing(pi.iip_y.v_g1, g2).0;
        if gt_B != gt_y {
            return false;
        }
    }
    // Group-level Mul gadget: enforce P(τ) = A(τ)·B(τ) - C(τ)
    {
        let g2 = <Bn254 as Pairing>::G2::generator();
        let gt_ab = <Bn254 as Pairing>::pairing(pi.a_tau_1, pi.b_tau_2).0;
        let gt_p  = <Bn254 as Pairing>::pairing(pi.p_tau_1, g2).0;
        let gt_c  = <Bn254 as Pairing>::pairing(pi.c_tau_1, g2).0;

        if gt_ab != gt_p * gt_c {
            return false;
        }
    }

    let shape = dg.linear_shape(crs);
    let coords = match build_lv_coords(crs, dg, pi) {
        Some(c) => c,
        None => return false,
    };

    recover_sb_via_linear_check(&shape, &coords.0)
}
