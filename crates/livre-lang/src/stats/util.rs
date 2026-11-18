//! Helpers shared by the Statistics and Machine Learning Toolbox modules: dimension arguments
//! (`dim`, `vecdim`, `'all'`), the groups of elements a reduction runs over, and argument checks.

use crate::error::{rt, LangError};
use crate::nd::{self, Dims};
use crate::value::*;

pub type Res<T> = Result<T, LangError>;

/// Which dimensions a reduction runs over.
#[derive(Debug, Clone, PartialEq)]
pub enum DimSel {
    /// The first non-singleton dimension.
    Default,
    All,
    /// 1-based dimensions.
    Dims(Vec<usize>),
}

const DIM_MESSAGE: &str = "Dimension argument must be a positive integer scalar, a vector of unique positive integers, or 'all'.";

/// A dimension argument: a positive integer, a vector of unique ones, `'all'`, or `[]`.
pub fn parse_dim(v: &Value) -> Res<DimSel> {
    if v.is_text() {
        let t = v.text()?.to_lowercase();
        if !t.is_empty() && "all".starts_with(&t) {
            return Ok(DimSel::All);
        }
        return rt(DIM_MESSAGE);
    }
    let Value::Num(m) = v else { return rt(DIM_MESSAGE) };
    if m.is_empty() {
        return Ok(DimSel::Default);
    }
    let mut dims = Vec::new();
    for &d in &m.data {
        if d < 1.0 || d.fract() != 0.0 || !d.is_finite() || dims.contains(&(d as usize)) {
            return rt(DIM_MESSAGE);
        }
        dims.push(d as usize);
    }
    Ok(DimSel::Dims(dims))
}

/// The elements a reduction combines: one list of source positions (in ascending order) per
/// element of the result.
pub struct Groups {
    pub dims: Dims,
    pub idx: Vec<Vec<usize>>,
}

pub fn groups(dims: &[usize], sel: &DimSel) -> Groups {
    let reduced: Vec<usize> = match sel {
        DimSel::Default => vec![dims.iter().position(|&d| d != 1).unwrap_or(0)],
        DimSel::All => (0..dims.len()).collect(),
        DimSel::Dims(v) => v.iter().map(|d| d - 1).filter(|&d| d < dims.len()).collect(),
    };
    let mut out_dims: Vec<usize> = dims.to_vec();
    for &d in &reduced {
        out_dims[d] = 1;
    }
    let n_out: usize = out_dims.iter().product();
    let mut idx = vec![Vec::new(); n_out];
    let total: usize = dims.iter().product();
    if reduced.is_empty() {
        // Every element is its own group.
        for (k, g) in idx.iter_mut().enumerate() {
            g.push(k);
        }
        return Groups { dims: nd::normalize(out_dims), idx };
    }
    let mut out_stride = vec![1usize; dims.len()];
    for k in 1..dims.len() {
        out_stride[k] = out_stride[k - 1] * out_dims[k - 1];
    }
    let mut coord = vec![0usize; dims.len()];
    for lin in 0..total {
        let mut o = 0;
        for k in 0..dims.len() {
            if !reduced.contains(&k) {
                o += coord[k] * out_stride[k];
            }
        }
        idx[o].push(lin);
        for k in 0..dims.len() {
            coord[k] += 1;
            if coord[k] < dims[k] {
                break;
            }
            coord[k] = 0;
        }
    }
    Groups { dims: nd::normalize(out_dims), idx }
}

/// The groups for a call that also needs the reduced size (`mu` and `sigma` of `zscore`).
pub fn apply(m: &Matrix, sel: &DimSel, mut f: impl FnMut(&[f64]) -> f64) -> Matrix {
    let g = groups(&m.dims, sel);
    let mut buf = Vec::new();
    let data = g
        .idx
        .iter()
        .map(|ix| {
            buf.clear();
            buf.extend(ix.iter().map(|&i| m.data[i]));
            f(&buf)
        })
        .collect();
    Matrix::from_dims(g.dims, data)
}

/// Marks a double result as single precision when the input was.
pub fn finish(m: Matrix, single: bool) -> Matrix {
    if single { m.cast_to(Kind::Single) } else { m }
}

pub fn is_single(m: &Matrix) -> bool {
    m.kind == Kind::Single
}

/// The numeric (double or single) data argument `i`, or MATLAB's "not enough input" error.
pub fn data_arg<'a>(a: &'a [Value], i: usize) -> Res<&'a Matrix> {
    match a.get(i) {
        Some(Value::Num(m)) => Ok(m),
        Some(other) => rt(format!("Undefined function for input arguments of type '{}'.", other.class_name())),
        None => rt("Not enough input arguments."),
    }
}

pub fn is_empty_00(m: &Matrix) -> bool {
    m.dims == [0, 0]
}

pub fn mean_of(v: &[f64]) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.iter().sum::<f64>() / v.len() as f64
}

/// The median of `v` (NaN when empty); `v` is reordered.
pub fn median_of(v: &mut Vec<f64>) -> f64 {
    if v.is_empty() || v.iter().any(|x| x.is_nan()) {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 { v[n / 2] } else { (v[n / 2 - 1] + v[n / 2]) / 2.0 }
}
