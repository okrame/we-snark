// src/gadgets/mod.rs

pub mod traits;
pub mod linear;
pub mod arithmetic;

pub use traits::{LVGadget, LVShapeBuilder};
pub use linear::{IIPGadget, NonZeroGadget, MaxDegGadget};
pub use arithmetic::MulGadget;
