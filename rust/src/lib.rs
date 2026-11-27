//! The search kernels of seedmatch: search rules 2 and 4 of `docs/design.md`.
//!
//! A code is a subject mask. Each table lookup reads `code & 15`, thus a code outside the masks
//! gives a wrong letter, but not a panic.

#![forbid(unsafe_code)]

pub mod hamming;
pub mod indel;
pub mod kernels;
pub mod letters;
#[cfg(feature = "python")]
mod python;

/// A slice has at most `isize::MAX` items, thus each length and each index fits in `i64`.
fn to_i64(value: usize) -> i64 {
    i64::try_from(value).unwrap_or(i64::MAX)
}

/// A slice has at most `isize::MAX` items, thus each length fits in `u64`.
fn to_u64(value: usize) -> u64 {
    u64::try_from(value).unwrap_or(u64::MAX)
}
