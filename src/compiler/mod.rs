// src/compiler/mod.rs
//! R1CS → QAP → LV “compiler” scaffolding.
//!
//! Goal for Phase 2/3:
//!   * define the data structures and function signatures we will need
//!     once we start consuming real `.r1cs` files;
//!   * keep the current Mul demo completely unchanged.

use ark_bn254::Fr;

use crate::gadgets::arithmetic::QAP;
use crate::scs::CRS;

/// Minimal shape information for an R1CS instance.
///
/// This is intentionally small; later we can enrich it with matrices A,B,C,
/// public-input indices, etc.
#[derive(Debug, Clone)]
pub struct R1CSShape {
    /// Number of constraints (rows of A,B,C).
    pub num_constraints: usize,
    /// Total number of variables (witness + public inputs).
    pub num_variables: usize,
    /// Number of *public* inputs.
    pub num_public_inputs: usize,
}

/// Sparse R1CS matrices in a compiler-friendly form.
///
/// We store A,B,C in *column-major* layout: for each variable j, the vector
/// `A_cols[j]` contains the (row_index, coefficient) pairs for which
/// A[row_index, j] != 0. This makes it easy to interpolate one polynomial
/// A_j(X) per variable when building a QAP.
#[derive(Debug, Clone)]
#[allow(non_snake_case)]
pub struct R1CSMatrices {
    pub num_constraints: usize,
    pub num_variables: usize,
    pub num_public_inputs: usize,
    pub A_cols: Vec<Vec<(usize, Fr)>>,
    pub B_cols: Vec<Vec<(usize, Fr)>>,
    pub C_cols: Vec<Vec<(usize, Fr)>>,
}

/// Placeholder loader: the *shape* and *matrix* representation we will
/// eventually populate from an actual `.r1cs` parser (ark-relations,
/// ark-circom, or a custom loader).
///
/// Keeping this function in place makes the rest of the Phase-3 API compile
/// even before the real parsing work is done.
pub fn load_r1cs(_path: &str) -> (R1CSShape, R1CSMatrices) {
    unimplemented!("load_r1cs: wire this to your chosen R1CS parser");
}


/// Result of the “algebraic” part of the compiler:
/// R1CS → QAP polynomials.
#[derive(Debug, Clone)]
pub struct CompiledQAP {
    pub shape: R1CSShape,
    pub qap: QAP,
}

impl CompiledQAP {
    /// Tiny helper for the current demo circuit w = [x,y,z,1].
    ///
    /// This is the “single gate” QAP. It just wraps `QAP::for_mul`
    /// so that the rest of the code can pretend it came from an R1CS.
    pub fn from_mul_demo(w: &[Fr]) -> Self {
        assert_eq!(
            w.len(),
            4,
            "mul demo expects witness [x,y,z,1] with 4 coordinates"
        );
        let qap = QAP::for_mul(w);
        let shape = R1CSShape {
            num_constraints: 1,
            num_variables: 3,    // x, y, z (ignoring the constant 1 slot)
            num_public_inputs: 0,
        };
        CompiledQAP { shape, qap }
    }

    /// Real Phase-3 entry point:
    ///
    /// * parses `path` into an `R1CSShape` + sparse matrices;
    /// * turns them into a QAP via `QAP::from_r1cs`.
    ///
    /// `crs` is currently unused here but kept in the signature so you
    /// can later add degree/size sanity checks if you want.
    pub fn from_r1cs(_crs: &CRS, path: &str) -> Self {
        let (shape, matrices) = load_r1cs(path);
        let qap = QAP::from_r1cs(&matrices);
        CompiledQAP { shape, qap }
    }
}


/// LV system-parameter helper.
///
/// This is where we centralise the choice of the MaxDeg bound `d`
/// for the witness polynomial B(X) and any other global variables that
/// depend on the circuit size.
#[derive(Debug, Clone, Copy)]
pub struct LVSystemParams {
    /// Degree bound for the B(X) polynomial in the IIP gadget.
    /// For the Mul demo we keep the old choice d = n − 1.
    pub d_bound: usize,
}

impl LVSystemParams {
    /// Old hard-coded choice used in the Mul demo:
    ///   * domain size n = 4,
    ///   * witness polynomial B(X) has degree ≤ n − 1.
    pub fn for_mul_demo(crs: &CRS) -> Self {
        LVSystemParams {
            d_bound: crs.n - 1,
        }
    }

    /// Derive LV system parameters from the R1CS shape.
    ///
    /// For now we take a very simple policy:
    ///   * degree bound d = m - 1 where m = num_constraints;
    ///   * require that the CRS supports degree ≥ d.
    ///
    /// This makes the IIP witness polynomial B(X) “fit inside” the CRS
    /// and is consistent with the original mul demo (where m = 1).
    pub fn from_r1cs(crs: &CRS, shape: &R1CSShape) -> Self {
        let m = shape.num_constraints;
        assert!(
            m > 0,
            "LVSystemParams::from_r1cs: zero-constraint shapes are not supported yet"
        );

        let d_bound = m.saturating_sub(1);

        assert!(
            d_bound <= crs.n - 1,
            "CRS too small: need degree bound d <= {}, but CRS.n = {}",
            d_bound,
            crs.n
        );

        LVSystemParams { d_bound }
    }
}
