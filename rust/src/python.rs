//! The module `seedmatch.search._kernels`. Each function releases the GIL for its full search.

use numpy::{PyArray1, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::kernels::{self, OffsetsError};
use crate::letters::{Compat, MASKS};

type Array<'py> = Bound<'py, PyArray1<i64>>;

fn compat_table(compat: &PyReadonlyArray2<'_, bool>) -> PyResult<Compat> {
    let view = compat.as_array();
    let (rows, columns) = view.dim();
    if (rows, columns) != (MASKS, MASKS) {
        return Err(PyValueError::new_err(format!(
            "compat has shape ({rows}, {columns}); expected ({MASKS}, {MASKS})"
        )));
    }
    let mut table = [[false; MASKS]; MASKS];
    for ((p, s), &value) in view.indexed_iter() {
        table[p][s] = value;
    }
    Ok(table)
}

fn value_error(error: OffsetsError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

#[pyfunction]
#[pyo3(signature = (codes, offsets, pattern, compat, max_mm, min_mm, /))]
fn hamming_counts<'py>(
    py: Python<'py>,
    codes: PyReadonlyArray1<'py, u8>,
    offsets: PyReadonlyArray1<'py, i64>,
    pattern: PyReadonlyArray1<'py, u8>,
    compat: PyReadonlyArray2<'py, bool>,
    max_mm: u64,
    min_mm: u64,
) -> PyResult<Array<'py>> {
    let compat = compat_table(&compat)?;
    let (codes, offsets, pattern) = (codes.as_slice()?, offsets.as_slice()?, pattern.as_slice()?);
    let counts = py
        .detach(|| kernels::hamming_counts(codes, offsets, pattern, &compat, max_mm, min_mm))
        .map_err(value_error)?;
    Ok(PyArray1::from_vec(py, counts))
}

#[pyfunction]
#[pyo3(signature = (codes, offsets, pattern, compat, max_mm, min_mm, /))]
fn hamming_ranges<'py>(
    py: Python<'py>,
    codes: PyReadonlyArray1<'py, u8>,
    offsets: PyReadonlyArray1<'py, i64>,
    pattern: PyReadonlyArray1<'py, u8>,
    compat: PyReadonlyArray2<'py, bool>,
    max_mm: u64,
    min_mm: u64,
) -> PyResult<(Array<'py>, Array<'py>)> {
    let compat = compat_table(&compat)?;
    let (codes, offsets, pattern) = (codes.as_slice()?, offsets.as_slice()?, pattern.as_slice()?);
    let (index, starts) = py
        .detach(|| kernels::hamming_ranges(codes, offsets, pattern, &compat, max_mm, min_mm))
        .map_err(value_error)?;
    Ok((
        PyArray1::from_vec(py, index),
        PyArray1::from_vec(py, starts),
    ))
}

#[pyfunction]
#[pyo3(signature = (codes, offsets, pattern, compat, max_mm, /))]
fn indel_counts<'py>(
    py: Python<'py>,
    codes: PyReadonlyArray1<'py, u8>,
    offsets: PyReadonlyArray1<'py, i64>,
    pattern: PyReadonlyArray1<'py, u8>,
    compat: PyReadonlyArray2<'py, bool>,
    max_mm: u64,
) -> PyResult<Array<'py>> {
    let compat = compat_table(&compat)?;
    let (codes, offsets, pattern) = (codes.as_slice()?, offsets.as_slice()?, pattern.as_slice()?);
    let counts = py
        .detach(|| kernels::indel_counts(codes, offsets, pattern, &compat, max_mm))
        .map_err(value_error)?;
    Ok(PyArray1::from_vec(py, counts))
}

#[pyfunction]
#[pyo3(signature = (codes, offsets, pattern, compat, max_mm, /))]
fn indel_ranges<'py>(
    py: Python<'py>,
    codes: PyReadonlyArray1<'py, u8>,
    offsets: PyReadonlyArray1<'py, i64>,
    pattern: PyReadonlyArray1<'py, u8>,
    compat: PyReadonlyArray2<'py, bool>,
    max_mm: u64,
) -> PyResult<(Array<'py>, Array<'py>, Array<'py>)> {
    let compat = compat_table(&compat)?;
    let (codes, offsets, pattern) = (codes.as_slice()?, offsets.as_slice()?, pattern.as_slice()?);
    let (index, starts, ends) = py
        .detach(|| kernels::indel_ranges(codes, offsets, pattern, &compat, max_mm))
        .map_err(value_error)?;
    Ok((
        PyArray1::from_vec(py, index),
        PyArray1::from_vec(py, starts),
        PyArray1::from_vec(py, ends),
    ))
}

#[pymodule]
fn _kernels(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(hamming_counts, module)?)?;
    module.add_function(wrap_pyfunction!(hamming_ranges, module)?)?;
    module.add_function(wrap_pyfunction!(indel_counts, module)?)?;
    module.add_function(wrap_pyfunction!(indel_ranges, module)?)?;
    Ok(())
}
