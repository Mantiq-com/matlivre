//! Statistics and Machine Learning Toolbox. Each area has its own module with a `NAMES` list and a
//! `call` dispatcher; `call` here routes a name to its module (see docs/stats-coverage-plan.md).

pub mod corr;
pub mod descr;
pub mod dist;
pub mod dist_matrix;
pub mod errtab;
pub mod hier;
pub mod kmeans;
pub mod htest;
pub mod pca;
pub mod regress;
pub mod sampling;
pub mod rnd;
pub mod util;

use crate::builtins::Ret;
use crate::interp::Interpreter;
use crate::value::*;

/// The function in `name`, if one of the stats modules provides it. Results of functions that do not
/// manage their class follow single precision input (the modules cast their own).
pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    if descr::NAMES.contains(&name) {
        return descr::call(it, name, a, nargout);
    }
    if dist::NAMES.contains(&name) {
        return dist::call(it, name, a, nargout);
    }
    if htest::NAMES.contains(&name) {
        return htest::call(it, name, a, nargout);
    }
    if dist_matrix::NAMES.contains(&name) {
        return dist_matrix::call(it, name, a, nargout);
    }
    if hier::NAMES.contains(&name) {
        return hier::call(it, name, a, nargout);
    }
    if pca::NAMES.contains(&name) {
        return pca::call(it, name, a, nargout);
    }
    if regress::NAMES.contains(&name) {
        return regress::call(it, name, a, nargout);
    }
    if kmeans::NAMES.contains(&name) {
        return kmeans::call(it, name, a, nargout);
    }
    if sampling::NAMES.contains(&name) {
        return sampling::call(it, name, a, nargout);
    }
    if rnd::NAMES.contains(&name) {
        return rnd::call(it, name, a, nargout);
    }
    if corr::NAMES.contains(&name) {
        return corr::call(it, name, a, nargout);
    }
    Ok(None)
}

/// Every name the stats modules provide.
pub fn names() -> Vec<&'static str> {
    [descr::NAMES, corr::NAMES, dist::NAMES, rnd::NAMES, htest::NAMES, dist_matrix::NAMES, hier::NAMES, pca::NAMES, regress::NAMES, kmeans::NAMES, sampling::NAMES].concat()
}
