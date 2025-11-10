use crate::iip::{IIPDigest, IIPProof, iip_verify};
use crate::nonzero::{NonZeroProof, nonzero_verify};
use crate::scs::CRS;
use ark_bn254::Bn254;
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ec::CurveGroup;

pub struct LVDigest {
    pub iip: IIPDigest,
    pub one_idx: usize,
}

pub struct LVProof {
    pub iip: IIPProof,
    pub nz: NonZeroProof,
}

/// Number of GT-coordinates we use in A_LV · π = b_LV.
pub const LV_NUM_COORDS: usize = 10;

/// A_LV and b_LV as described above.
/// - a[i][j] ∈ {-1,0,1} describes exponent α_{i,j} on coordinate c_j in equation i.
/// - b[i] ∈ GT is the RHS constant for equation i.
pub struct LVShape {
    pub rows: usize,
    pub a: [[i8; LV_NUM_COORDS]; 4], // here rows=4 fixed, you can generalize later
    pub b: [PairingOutput<Bn254>; 4],
}

impl LVDigest {
    /// Build the A_LV and b_LV used by the LV verifier, in the abstract
    /// GT-coordinate basis c_0..c_9 described in the documentation.
    pub fn linear_shape(&self, crs: &CRS) -> LVShape {
        use ark_bn254::Fr;
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
        let gt_one = PairingOutput::<Bn254>::one();
        let gt_const = <Bn254 as Pairing>::pairing(
            <Bn254 as Pairing>::G1::generator(),
            <Bn254 as Pairing>::G2::generator(),
        );

        let b = [gt_one, gt_one, gt_one, gt_const];

        LVShape { rows, a, b }
    }
}


pub fn lv_verify(_crs: &CRS, dg: &LVDigest, pi: &LVProof) -> bool {
    // 1) IIP checks
    if !iip_verify(&dg.iip, &pi.iip) { return false; }
    // 2) NonZero (actually “slot equals 1”) check
    if !nonzero_verify(&_crs, &pi.nz, dg.one_idx) { return false; }
    true
}
