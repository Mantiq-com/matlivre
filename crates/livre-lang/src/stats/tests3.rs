//! `jbtest` and `dwtest`.

use super::util::*;
use crate::builtins::Ret;
use crate::error::rt;
use crate::interp::Interpreter;
use crate::value::*;

pub const NAMES: &[&str] = &["jbtest", "dwtest", "chi2gof"];

const JB_ALPHAS: [f64; 17] = [0.001, 0.0016681, 0.0027826, 0.0046416, 0.0077426, 0.01, 0.012915, 0.021544, 0.025, 0.035938, 0.05, 0.059948, 0.1, 0.15, 0.2, 0.25, 0.50];
const JB_SIZES: [f64; 32] = [
    3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0, 300.0, 400.0, 500.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0,
];
const JB_TABLE: [[f64; 17]; 32] = [
    [0.5312, 0.5312, 0.5312, 0.5312, 0.5312, 0.5312, 0.5311, 0.531, 0.5309, 0.5305, 0.5297, 0.529, 0.5251, 0.5176, 0.5074, 0.4946, 0.4063],
    [0.9606, 0.959, 0.9564, 0.952, 0.9448, 0.9396, 0.9329, 0.9133, 0.9056, 0.8817, 0.8519, 0.8316, 0.7553, 0.6721, 0.6303, 0.5947, 0.4739],
    [1.8289, 1.8053, 1.7727, 1.7276, 1.6661, 1.6275, 1.5828, 1.4723, 1.4343, 1.3294, 1.2185, 1.1516, 0.9442, 0.7945, 0.7302, 0.6878, 0.5285],
    [10.9719, 9.843, 8.6751, 7.4815, 6.2927, 5.7077, 5.135, 4.0456, 3.7474, 3.0691, 2.5239, 2.2555, 1.6231, 1.2821, 1.1235, 1.0198, 0.6951],
    [19.5425, 16.7207, 14.0454, 11.5527, 9.2814, 8.2365, 7.2613, 5.5238, 5.0729, 4.0773, 3.2985, 2.9215, 2.0533, 1.5965, 1.3779, 1.2336, 0.7916],
    [25.0722, 21.0001, 17.2981, 13.9762, 11.0602, 9.7531, 8.5489, 6.4401, 5.8998, 4.7176, 3.8011, 3.3596, 2.3505, 1.8239, 1.5631, 1.3885, 0.8577],
    [28.4885, 23.6332, 19.2818, 15.4635, 12.1658, 10.7058, 9.3658, 7.0389, 6.4463, 5.1501, 4.1494, 3.6684, 2.5707, 1.9986, 1.7063, 1.5079, 0.9071],
    [30.6106, 25.2501, 20.5242, 16.4088, 12.8805, 11.3263, 9.9063, 7.4472, 6.8222, 5.4578, 4.4039, 3.8968, 2.7431, 2.139, 1.8216, 1.604, 0.9457],
    [31.9343, 26.2695, 21.3105, 17.0176, 13.3572, 11.7488, 10.2818, 7.7395, 7.0942, 5.6856, 4.5973, 4.0734, 2.8827, 2.2547, 1.9172, 1.6835, 0.9771],
    [32.7514, 26.9025, 21.8016, 17.4145, 13.6749, 12.0355, 10.5399, 7.952, 7.2949, 5.8598, 4.7481, 4.2126, 2.9987, 2.3524, 1.998, 1.7508, 1.0033],
    [33.2273, 27.2686, 22.103, 17.6594, 13.8848, 12.2302, 10.7201, 8.1098, 7.4469, 5.9945, 4.8689, 4.3258, 3.0973, 2.4361, 2.0674, 1.8085, 1.0256],
    [33.4825, 27.4765, 22.2863, 17.8172, 14.0345, 12.3739, 10.8586, 8.234, 7.5661, 6.1045, 4.9697, 4.4214, 3.1834, 2.5097, 2.1283, 1.8592, 1.0449],
    [33.5009, 27.5458, 22.379, 17.9431, 14.1788, 12.5255, 11.0182, 8.3968, 7.7287, 6.262, 5.1203, 4.5688, 3.3246, 2.6316, 2.2297, 1.9434, 1.0768],
    [33.261, 27.3738, 22.2848, 17.9209, 14.2092, 12.5801, 11.0884, 8.4921, 7.8286, 6.3696, 5.2305, 4.6803, 3.4374, 2.7298, 2.3115, 2.0114, 1.1023],
    [32.8742, 27.0975, 22.1135, 17.8305, 14.1911, 12.5841, 11.1155, 8.552, 7.8938, 6.4462, 5.3134, 4.7662, 3.5292, 2.8104, 2.3789, 2.0674, 1.1231],
    [32.3418, 26.7278, 21.877, 17.6977, 14.1283, 12.5542, 11.1112, 8.585, 7.9361, 6.5028, 5.3796, 4.8378, 3.6071, 2.8788, 2.4361, 2.1151, 1.1408],
    [31.8884, 26.399, 21.6407, 17.551, 14.0491, 12.5067, 11.0882, 8.6004, 7.9589, 6.5429, 5.4314, 4.8957, 3.673, 2.9372, 2.4851, 2.1557, 1.1557],
    [30.6089, 25.4748, 21.01, 17.146, 13.8271, 12.3565, 10.9976, 8.602, 7.9816, 6.6071, 5.5277, 5.0071, 3.8025, 3.0523, 2.5819, 2.2363, 1.1852],
    [29.429, 24.6043, 20.3979, 16.7425, 13.5878, 12.1804, 10.8779, 8.5719, 7.9718, 6.6387, 5.5919, 5.0863, 3.8987, 3.1384, 2.6543, 2.2964, 1.2072],
    [28.4225, 23.8582, 19.8658, 16.3846, 13.3664, 12.0139, 10.7589, 8.5307, 7.9496, 6.6571, 5.6408, 5.1474, 3.9732, 3.205, 2.7106, 2.3436, 1.2243],
    [27.514, 23.1914, 19.3922, 16.0635, 13.1643, 11.8602, 10.6494, 8.4904, 7.9259, 6.6688, 5.6783, 5.1957, 4.0327, 3.2584, 2.7559, 2.3812, 1.2382],
    [26.0239, 22.0764, 18.5877, 15.5123, 12.8124, 11.5923, 10.4542, 8.4114, 7.8751, 6.6803, 5.7343, 5.2681, 4.1224, 3.3402, 2.825, 2.4388, 1.2591],
    [24.8353, 21.1758, 17.9384, 15.0633, 12.5206, 11.3653, 10.2839, 8.3384, 7.8276, 6.6842, 5.7728, 5.3194, 4.1873, 3.3988, 2.8749, 2.4808, 1.2744],
    [23.0771, 19.8542, 16.9669, 14.3843, 12.0743, 11.0157, 10.0226, 8.2266, 7.7526, 6.6877, 5.8248, 5.3889, 4.2748, 3.479, 2.9434, 2.5382, 1.2956],
    [21.817, 18.9068, 16.2666, 13.8903, 11.7469, 10.7604, 9.8318, 8.1454, 7.6984, 6.6882, 5.8581, 5.4336, 4.332, 3.5318, 2.9886, 2.576, 1.3097],
    [19.5539, 17.1748, 14.9865, 12.9765, 11.1441, 10.2914, 9.4852, 8.0027, 7.604, 6.6858, 5.9096, 5.5043, 4.4246, 3.618, 3.0631, 2.6388, 1.3334],
    [18.6707, 16.4879, 14.4745, 12.6143, 10.9072, 10.1106, 9.3521, 7.9501, 7.5691, 6.6847, 5.9282, 5.5296, 4.458, 3.65, 3.0907, 2.6624, 1.3421],
    [18.0126, 15.983, 14.099, 12.3493, 10.7362, 9.9803, 9.2587, 7.9127, 7.5441, 6.6832, 5.9408, 5.5471, 4.4814, 3.6718, 3.1099, 2.6787, 1.3484],
    [17.519, 15.6036, 13.8149, 12.1501, 10.6072, 9.8811, 9.1869, 7.8828, 7.5239, 6.68, 5.9482, 5.5587, 4.4979, 3.6876, 3.1237, 2.6904, 1.353],
    [17.1305, 15.3029, 13.5883, 11.9883, 10.5059, 9.8072, 9.1333, 7.8611, 7.5089, 6.6789, 5.9546, 5.5677, 4.5105, 3.6997, 3.1345, 2.6996, 1.3564],
    [16.8158, 15.06, 13.4102, 11.8682, 10.4296, 9.7478, 9.0911, 7.8446, 7.4979, 6.6777, 5.96, 5.5755, 4.5214, 3.7101, 3.1433, 2.7072, 1.3595],
    [16.5585, 14.8633, 13.2657, 11.7666, 10.3636, 9.6981, 9.0549, 7.8296, 7.4872, 6.6763, 5.9635, 5.5806, 4.5289, 3.7173, 3.1499, 2.7129, 1.3619],
];

/// Shape-preserving cubic interpolation (the slopes of Fritsch and Carlson that MATLAB's `pchip`
/// uses) of `(x, y)` at `q`; `x` increases strictly.
fn pchip_at(x: &[f64], y: &[f64], q: f64) -> f64 {
    let n = x.len();
    let h: Vec<f64> = (0..n - 1).map(|i| x[i + 1] - x[i]).collect();
    let del: Vec<f64> = (0..n - 1).map(|i| (y[i + 1] - y[i]) / h[i]).collect();
    let mut d = vec![0.0; n];
    if n == 2 {
        d[0] = del[0];
        d[1] = del[0];
    } else {
        for k in 1..n - 1 {
            if del[k - 1] * del[k] > 0.0 {
                let hs = h[k - 1] + h[k];
                let w1 = (h[k - 1] + hs) / (3.0 * hs);
                let w2 = (hs + h[k]) / (3.0 * hs);
                let (dmax, dmin) = (del[k - 1].abs().max(del[k].abs()), del[k - 1].abs().min(del[k].abs()));
                d[k] = dmin / (w1 * (del[k - 1] / dmax) + w2 * (del[k] / dmax));
            }
        }
        let end = |h0: f64, h1: f64, d0: f64, d1: f64| {
            let mut s = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
            if s.signum() != d0.signum() {
                s = 0.0;
            } else if d0.signum() != d1.signum() && s.abs() > 3.0 * d0.abs() {
                s = 3.0 * d0;
            }
            s
        };
        d[0] = end(h[0], h[1], del[0], del[1]);
        d[n - 1] = end(h[n - 2], h[n - 3], del[n - 2], del[n - 3]);
    }
    let mut k = 0;
    while k + 2 < n && q >= x[k + 1] {
        k += 1;
    }
    let t = q - x[k];
    let c = (3.0 * del[k] - 2.0 * d[k] - d[k + 1]) / h[k];
    let b = (d[k] - 2.0 * del[k] + d[k + 1]) / (h[k] * h[k]);
    y[k] + t * (d[k] + t * (c + t * b))
}

/// The critical values of the Jarque-Bera statistic for sample size `n` at the tabulated levels.
fn jb_critical_values(n: f64) -> Vec<f64> {
    // Sizes enter as 1/n; the last row is the limit (chi-square with two degrees of freedom).
    let mut inv: Vec<f64> = vec![0.0];
    inv.extend(JB_SIZES.iter().rev().map(|s| 1.0 / s));
    (0..17)
        .map(|j| {
            let mut col = vec![((-2.0 * JB_ALPHAS[j].ln()) * 1e4).round() / 1e4];
            col.extend(JB_TABLE.iter().rev().map(|row| row[j]));
            pchip_at(&inv, &col, 1.0 / n)
        })
        .collect()
}

fn optional<'a>(a: &'a [Value], i: usize) -> Option<&'a Value> {
    match a.get(i) {
        Some(Value::Num(m)) if m.is_empty() => None,
        other => other,
    }
}

fn jbtest(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let x = match a.first() {
        Some(Value::Num(m)) => m,
        Some(_) => return rt("Input X must be a vector of real values."),
        None => return rt("Not enough input arguments."),
    };
    let is_vec = x.dims.len() == 2 && (x.dims[0] == 1 || x.dims[1] == 1);
    if !is_vec {
        return rt("Input X must be a vector of real values.");
    }
    let xs: Vec<f64> = x.data.iter().copied().filter(|v| !v.is_nan()).collect();
    let n = xs.len();
    if n < 2 {
        return rt("Sample vector X must have at least 2 valid observations.");
    }
    if n == 2 {
        return Ok(Some(vec![Value::num(0.0), Value::num(1.0), Value::num(1.0 / 3.0), Value::num(f64::NAN)]));
    }
    let alpha = match optional(a, 1) {
        Some(Value::Num(m)) if m.data.len() == 1 && 0.0 < m.data[0] && m.data[0] < 1.0 => m.data[0],
        Some(_) => return rt("Significance level ALPHA must be a scalar between 0 and 1."),
        None => 0.05,
    };
    let mctol = match optional(a, 2) {
        Some(Value::Num(m)) if m.data.len() == 1 && m.data[0] > 0.0 => Some(m.data[0]),
        Some(_) => return rt("Monte-Carlo standard error tolerance MCTOL must be a positive scalar value."),
        None => None,
    };
    let nf = n as f64;
    let mean = xs.iter().sum::<f64>() / nf;
    let sd1 = (xs.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / nf).sqrt();
    let z: Vec<f64> = xs.iter().map(|v| (v - mean) / sd1).collect();
    let skew = z.iter().map(|v| v.powi(3)).sum::<f64>() / nf;
    let kurt = z.iter().map(|v| v.powi(4)).sum::<f64>() / nf - 3.0;
    let jb = nf * (skew * skew / 6.0 + kurt * kurt / 24.0);
    let (cv, p);
    if let Some(tol) = mctol {
        let (c, pp) = jb_monte_carlo(it, jb, n, alpha, tol)?;
        cv = c;
        p = pp;
    } else {
        if alpha < JB_ALPHAS[0] || JB_ALPHAS[16] < alpha {
            return rt("Significance level ALPHA is outside the range of the tabulated values. Use a value in the interval [0.001, 0.5], or use the MCTOL input argument.");
        }
        let cvs = jb_critical_values(nf);
        cv = pchip_at(&JB_ALPHAS, &cvs, alpha);
        let mut pv = f64::NAN;
        if nargout > 1 {
            if jb < cvs[16] {
                it.emit_warning("P is greater than the largest tabulated value, returning 0.5.".into(), "stats:jbtest:PTooBig".into());
                pv = JB_ALPHAS[16];
            } else if cvs[0] < jb {
                it.emit_warning("P is less than the smallest tabulated value, returning 0.001.".into(), "stats:jbtest:PTooSmall".into());
                pv = JB_ALPHAS[0];
            } else if jb == cvs[0] {
                pv = JB_ALPHAS[0];
            } else if jb.is_nan() {
                pv = 0.0;
            } else {
                let i = cvs.iter().position(|&c| jb > c).unwrap();
                let (mut lo, mut hi) = (JB_ALPHAS[i - 1], JB_ALPHAS[i]);
                let f = |t: f64| pchip_at(&JB_ALPHAS, &cvs, t) - jb;
                let flo = f(lo);
                for _ in 0..200 {
                    let mid = 0.5 * (lo + hi);
                    if (f(mid) > 0.0) == (flo > 0.0) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                pv = 0.5 * (lo + hi);
            }
        }
        p = pv;
    }
    let h = if jb > cv { 1.0 } else { 0.0 } + if jb.is_nan() { 1.0 } else { 0.0 };
    let single = is_single(x);
    let sv = |v: f64| Value::Num(finish(Matrix::scalar(v), single));
    Ok(Some(vec![Value::num(h), sv(p), sv(jb), Value::num(cv)]))
}

fn jb_monte_carlo(it: &mut Interpreter, jb: f64, n: usize, alpha: f64, tol: f64) -> Res<(f64, f64)> {
    let vartol = tol * tol;
    let (mut crit, mut p, mut reps_tot) = (0.0, 0.0, 0.0f64);
    let mut reps_min = 1000.0f64;
    loop {
        let old = reps_tot;
        let reps = (reps_min - old).ceil().max(1.0) as usize;
        let mut stats = Vec::with_capacity(reps);
        for _ in 0..reps {
            let r = it.call_function("randn", vec![Value::num(1.0), Value::num(n as f64)], 1)?;
            let Value::Num(m) = &r[0] else { return rt("randn") };
            let nf = n as f64;
            let mean = m.data.iter().sum::<f64>() / nf;
            let s2 = m.data.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / nf;
            let skew = m.data.iter().map(|v| (v - mean).powi(3)).sum::<f64>() / (s2.powf(1.5) * nf);
            let kurt = m.data.iter().map(|v| (v - mean).powi(4)).sum::<f64>() / (s2 * s2 * nf) - 3.0;
            stats.push(nf * (skew * skew / 6.0 + kurt * kurt / 24.0));
        }
        let q = it.call_function("prctile", vec![Value::Num(Matrix::from_dims(vec![reps, 1], stats.clone())), Value::num(100.0 * (1.0 - alpha))], 1)?;
        let crit_mc = q[0].scalar()?;
        let p_mc = stats.iter().filter(|&&s| s > jb).count() as f64 / reps as f64;
        reps_tot = old + reps as f64;
        crit = (old * crit + reps as f64 * crit_mc) / reps_tot;
        p = (old * p + reps as f64 * p_mc) / reps_tot;
        let sepsq = (p * (1.0 - p) / reps_tot).max(1.0 / (reps_tot * reps_tot));
        if sepsq < vartol {
            break;
        }
        reps_min = 1.2 * (reps_tot * sepsq) / vartol;
    }
    if jb.is_nan() {
        p = 0.0;
    }
    Ok((crit, p))
}

// ---------------------------------------------------------------------------------------------
// dwtest

type Mat = Vec<Vec<f64>>;

fn matmul(a: &Mat, b: &Mat) -> Mat {
    let (n, k, m) = (a.len(), b.len(), b[0].len());
    let mut c = vec![vec![0.0; m]; n];
    for i in 0..n {
        for l in 0..k {
            let v = a[i][l];
            for j in 0..m {
                c[i][j] += v * b[l][j];
            }
        }
    }
    c
}

fn transpose(a: &Mat) -> Mat {
    (0..a[0].len()).map(|j| a.iter().map(|r| r[j]).collect()).collect()
}

fn trace(a: &Mat) -> f64 {
    (0..a.len()).map(|i| a[i][i]).sum()
}

fn inverse(a: &Mat) -> Mat {
    let n = a.len();
    let mut m: Mat = (0..n).map(|i| a[i].iter().copied().chain((0..n).map(|j| if i == j { 1.0 } else { 0.0 })).collect()).collect();
    for c in 0..n {
        let piv = (c..n).max_by(|&i, &j| m[i][c].abs().partial_cmp(&m[j][c].abs()).unwrap()).unwrap();
        m.swap(c, piv);
        let d = m[c][c];
        for v in m[c].iter_mut() {
            *v /= d;
        }
        for r in 0..n {
            if r != c {
                let f = m[r][c];
                if f != 0.0 {
                    for j in 0..2 * n {
                        let t = m[c][j];
                        m[r][j] -= f * t;
                    }
                }
            }
        }
    }
    m.into_iter().map(|r| r[n..].to_vec()).collect()
}

/// P(DW < dw) from the eigenvalues `lam` (ascending) of the quadratic form, by numerical
/// integration of the characteristic function (Pan's method as described by Farebrother).
fn pan_cdf(lam: &[f64], dw: f64, n: usize) -> f64 {
    let m = lam.len() as i64;
    let Some(first) = lam.iter().position(|&l| l >= dw) else { return 0.0 };
    let mut mu = first as i64;
    let mut h = m - mu;
    let (d, mut j1, j2, mut j3, mut j4, k);
    if h > mu {
        d = 2;
        h = mu;
        k = -1.0;
        j1 = 0;
        j2 = 2;
        j3 = 3;
        j4 = 1;
    } else {
        d = -2;
        mu += 1;
        j1 = m - 2;
        j2 = m - 1;
        j3 = m + 1;
        j4 = m;
        k = 1.0;
    }
    let step = std::f64::consts::FRAC_PI_2 / n as f64;
    let mut sum = (k + 1.0) / 2.0;
    let mut sgn = k / n as f64;
    let at = |i: i64| lam[(i - 1) as usize];
    let mut start = j2;
    let mut l1 = h - 2 * (h / 2);
    while l1 >= 0 {
        let mut l2 = start;
        while (d > 0 && l2 <= mu) || (d < 0 && l2 >= mu) {
            let top = at(j4);
            let other = if l2 == 0 { dw } else { at(l2) };
            let u = 0.5 * (top + other);
            let v = 0.5 * (top - other);
            let mut acc = 0.0;
            let mut i = 1;
            while i <= 2 * n - 1 {
                let y = u - v * (i as f64 * step).cos();
                let mut prod = 1.0;
                for t in (1..=j1).chain(j3..=m) {
                    prod *= (y - dw) / (y - at(t));
                }
                acc += prod.abs().sqrt();
                i += 2;
            }
            sgn = -sgn;
            sum += sgn * acc;
            j1 += d;
            j3 += d;
            j4 += d;
            l2 += d;
        }
        if d == 2 {
            j3 -= 1;
        } else {
            j1 += 1;
        }
        start = 0;
        mu = 0;
        l1 -= 1;
    }
    sum
}

fn dw_approx(it: &mut Interpreter, x: &Mat, dw: f64) -> Res<f64> {
    let (n, p) = (x.len(), x[0].len());
    let xt = transpose(x);
    let q1 = inverse(&matmul(&xt, x));
    let mut b: Mat = Vec::with_capacity(n);
    b.push((0..p).map(|j| x[0][j] - x[1][j]).collect());
    for r in 1..n - 1 {
        b.push((0..p).map(|j| 2.0 * x[r][j] - x[r - 1][j] - x[r + 1][j]).collect());
    }
    b.push((0..p).map(|j| -x[n - 2][j] + x[n - 1][j]).collect());
    let c = matmul(&matmul(&xt, &b), &q1);
    let nu1 = 2.0 * (n as f64 - 1.0) - trace(&c);
    let btb = matmul(&matmul(&transpose(&b), &b), &q1);
    let nu2 = 2.0 * (3.0 * n as f64 - 4.0) - 2.0 * trace(&btb) + trace(&matmul(&c, &c));
    let df = n as f64 - p as f64;
    let mu = nu1 / df;
    let sigma = (2.0 / (df * (df + 2.0)) * (nu2 - nu1 * mu)).sqrt();
    let r = it.call_function("normcdf", vec![Value::num(dw), Value::num(mu), Value::num(sigma)], 1)?;
    r[0].scalar()
}

fn pick(given: &str, choices: &[&str], what: &str) -> Res<String> {
    let g = given.to_lowercase();
    if !g.is_empty() {
        let hits: Vec<&&str> = choices.iter().filter(|c| c.starts_with(&g)).collect();
        if hits.len() == 1 {
            return Ok(hits[0].to_string());
        }
    }
    let list: Vec<String> = choices.iter().map(|c| format!("'{c}'")).collect();
    rt(format!("'{given}' is not a valid value for the '{what}' argument.  Valid values are: {}.", list.join(", ")))
}

fn text_of(v: &Value) -> String {
    if v.is_text() {
        v.text().unwrap_or_default()
    } else {
        String::new()
    }
}

fn dwtest(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let _ = nargout;
    let (Some(Value::Num(r)), Some(Value::Num(xm))) = (a.first(), a.get(1)) else { return rt("Not enough input arguments.") };
    let (n, p) = (xm.dims[0], xm.dims[1]);
    let mut option = String::new();
    let mut alt = String::new();
    if a.len() > 2 {
        let v = &a[2];
        if !v.is_text() {
            return rt(format!("Input A of class {} and input B of class cell must be cell arrays of character vectors, unless one is a character vector.", v.class_name()));
        }
        let low = text_of(v).to_lowercase();
        if low == "exact" || low == "approximate" {
            option = low;
            if a.len() > 3 {
                alt = text_of(&a[3]);
            }
        } else if a.len() == 3 {
            return rt("The METHOD argument must be 'approximate' or 'exact'.");
        } else {
            let mut i = 2;
            while i < a.len() {
                let key = text_of(&a[i]).to_lowercase();
                let val = a.get(i + 1).map(text_of).unwrap_or_default();
                if !key.is_empty() && "method".starts_with(&key) {
                    option = val;
                } else if !key.is_empty() && "tail".starts_with(&key) {
                    alt = val;
                } else {
                    return rt(format!("Invalid parameter name: {}.", text_of(&a[i])));
                }
                i += 2;
            }
        }
    }
    let option = if option.is_empty() {
        if n < 400 { "exact".to_string() } else { "approximate".to_string() }
    } else {
        pick(&option, &["exact", "approximate"], "method")?
    };
    let alt = if alt.is_empty() { "both".to_string() } else { pick(&alt, &["left", "both", "right"], "tail")? };
    let d: Vec<f64> = r.data.windows(2).map(|w| w[1] - w[0]).collect();
    let dw = d.iter().map(|v| v * v).sum::<f64>() / r.data.iter().map(|v| v * v).sum::<f64>();
    let xrows: Mat = (0..n).map(|i| (0..p).map(|j| xm.data[j * n + i]).collect()).collect();
    let mut pdw = f64::NAN;
    let mut use_exact = option == "exact";
    if use_exact {
        let q = it.call_function("qr", vec![Value::Num(xm.clone())], 2)?;
        let Value::Num(qm) = &q[0] else { return rt("qr") };
        // G = F * Q(:, p+1:n), F the (n-1)-by-n difference matrix.
        let cols = n - p;
        let mut g = vec![0.0; (n - 1) * cols];
        for c in 0..cols {
            for i in 0..n - 1 {
                g[c * (n - 1) + i] = qm.data[(p + c) * n + i + 1] - qm.data[(p + c) * n + i];
            }
        }
        let sv = it.call_function("svd", vec![Value::Num(Matrix::from_dims(vec![n - 1, cols], g))], 1)?;
        let Value::Num(s) = &sv[0] else { return rt("svd") };
        let mut e: Vec<f64> = s.data.iter().map(|v| v * v).collect();
        e.sort_by(|a, b| a.partial_cmp(b).unwrap());
        pdw = pan_cdf(&e, dw, n);
        if pdw > 1.0 || pdw < 0.0 {
            it.emit_warning("Exact p-value is unavailable, using approximate p-value.".into(), "stats:pvaluedw:ExactUnavailable".into());
            use_exact = false;
        }
    }
    if !use_exact {
        pdw = dw_approx(it, &xrows, dw)?;
    }
    let pv = match alt.as_str() {
        "both" => 2.0 * pdw.min(1.0 - pdw),
        "left" => 1.0 - pdw,
        _ => pdw,
    };
    Ok(Some(vec![Value::num(pv), Value::num(dw)]))
}

pub fn call(it: &mut Interpreter, name: &str, a: &[Value], nargout: usize) -> Ret {
    match name {
        "jbtest" => jbtest(it, a, nargout),
        "chi2gof" => chi2gof(it, a, nargout),
        _ => dwtest(it, a, nargout),
    }
}

// ---------------------------------------------------------------------------------------------
// chi2gof

fn is_empty_val(v: &Value) -> bool {
    matches!(v, Value::Num(m) if m.is_empty())
}

fn row(v: Vec<f64>) -> Value {
    let n = v.len();
    Value::Num(Matrix::from_dims(vec![1, n], v))
}

fn is_vector_dims(d: &[usize]) -> bool {
    d.len() == 2 && (d[0] == 1 || d[1] == 1)
}

/// The bin counts and edges of `x` for a bin specification (nbins, centers or edges).
fn get_bins(it: &mut Interpreter, x: &[f64], freq: &[f64], kind: &str, spec: &Matrix) -> Res<(Vec<f64>, Vec<f64>)> {
    let lo = x.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut edges: Vec<f64>;
    let use_centers = kind != "edges";
    match kind {
        "nbins" => {
            let nb = spec.data[0];
            let (mut lo, mut hi) = if x.is_empty() { (0.0, 1.0) } else { (lo, hi) };
            if lo == hi {
                lo = lo - (nb / 2.0).floor() - 0.5;
                hi = hi + (nb / 2.0).ceil() - 0.5;
            }
            let bw = (hi - lo) / nb;
            let count = if nb < 0.0 { 0 } else { nb.floor() as usize + 1 };
            edges = (0..count).map(|k| lo + bw * k as f64).collect();
            if edges.is_empty() {
                return rt("Array indices must be positive integers or logical values.");
            }
            if count == 1 {
                return rt("First input IND and third input SZ must satisfy ALL(MAX(IND)<=SZ).");
            }
            let last = edges.len() - 1;
            edges[last] = hi;
        }
        "ctrs" => {
            let c = &spec.data;
            let mut bw: Vec<f64> = c.windows(2).map(|w| w[1] - w[0]).collect();
            if bw.iter().any(|&b| b <= 0.0) {
                return rt("The vector C of bin centers must be increasing.");
            }
            let l = *bw.last().unwrap_or(&f64::NAN);
            bw.push(l);
            edges = vec![c[0] - bw[0] / 2.0];
            edges.extend(c.iter().zip(&bw).map(|(ci, b)| ci + b / 2.0));
            if !x.is_empty() {
                edges[0] = edges[0].min(lo);
                let last = edges.len() - 1;
                edges[last] = edges[last].max(hi);
            }
        }
        _ => {
            edges = spec.data.clone();
            if edges.windows(2).any(|w| w[1] - w[0] <= 0.0) {
                return rt("The vector E of bin edges must be increasing.");
            }
            if !x.is_empty() && (lo < edges[0] || edges[edges.len() - 1] < hi) {
                it.emit_warning("The vector E of bin edges does not span the range of X.".into(), "stats:ecdfhist:DataNotSpanned".into());
            }
        }
    }
    let nb = edges.len() - 1;
    let mut counts = vec![0.0; nb];
    if use_centers {
        for e in edges.iter_mut() {
            *e += eps_of(*e);
        }
        for (v, w) in x.iter().zip(freq) {
            let k = edges[1..nb].iter().filter(|&&e| e <= *v).count();
            counts[k] += w;
        }
    } else {
        for (v, w) in x.iter().zip(freq) {
            if *v < edges[0] || *v > edges[nb] {
                continue;
            }
            let mut k = edges.iter().filter(|&&e| e <= *v).count().saturating_sub(1);
            if k >= nb {
                k = nb - 1;
            }
            counts[k] += w;
        }
    }
    Ok((counts, edges))
}

fn chi2gof(it: &mut Interpreter, a: &[Value], nargout: usize) -> Ret {
    let Some(Value::Num(xm)) = a.first() else { return rt("Not enough input arguments.") };
    if !is_vector_dims(&xm.dims) {
        return rt("X must be a vector of real values.");
    }
    let names = ["nbins", "ctrs", "edges", "cdf", "expected", "nparams", "emin", "frequency", "alpha"];
    let mut vals: Vec<Option<Value>> = vec![None; names.len()];
    if (a.len() - 1) % 2 != 0 {
        return rt("Wrong number of arguments.");
    }
    for pair in a[1..].chunks(2) {
        let key = text_of(&pair[0]).to_lowercase();
        let hits: Vec<usize> = names.iter().enumerate().filter(|(_, n)| !key.is_empty() && n.starts_with(&key)).map(|(i, _)| i).collect();
        let exact = names.iter().position(|n| *n == key);
        let idx = match (exact, hits.len()) {
            (Some(i), _) => i,
            (None, 1) => hits[0],
            (None, n) if n > 1 => return rt(format!("Ambiguous parameter name: {}.", text_of(&pair[0]))),
            _ => return rt(format!("Invalid parameter name: {}.", text_of(&pair[0]))),
        };
        vals[idx] = Some(pair[1].clone());
    }
    let get = |i: usize| -> Option<&Value> { vals[i].as_ref().filter(|v| !is_empty_val(v)) };
    let (nbins_v, ctrs_v, edges_v, cdf_v, exp_v, np_v, emin_v, freq_v, alpha_v) = (get(0), get(1), get(2), get(3), get(4), get(5), vals[6].as_ref(), get(7), vals[8].as_ref());
    if [nbins_v, ctrs_v, edges_v].iter().filter(|v| v.is_some()).count() > 1 {
        return rt("You cannot specify more than one of 'nbins', 'ctrs', and 'edges'.");
    }
    if cdf_v.is_some() && exp_v.is_some() {
        return rt("You cannot specify both 'cdf' and 'expected'.");
    }
    let freq_in: Option<Matrix> = match freq_v {
        Some(Value::Num(m)) => Some(m.clone()),
        _ => None,
    };
    if let Some(f) = &freq_in {
        if !is_vector_dims(&f.dims) || f.data.len() != xm.data.len() {
            return rt("The 'frequency' vector must have the same size as X.");
        }
        if f.data.iter().any(|&v| v < 0.0) {
            return rt("The 'frequency' values must be zero or positive.");
        }
    }
    let emin = match emin_v {
        None => 5.0,
        Some(Value::Num(m)) if m.data.len() == 1 && m.data[0] >= 0.0 && m.data[0] == m.data[0].round() => m.data[0],
        _ => return rt("The 'emin' value must be a non-negative integer."),
    };
    let mut nparams: Option<f64> = match np_v {
        None => None,
        Some(Value::Num(m)) if m.data.len() == 1 && m.data[0] >= 0.0 && m.data[0] == m.data[0].round() => Some(m.data[0]),
        _ => return rt("The 'nparams' value must be a non-negative integer."),
    };
    let alpha = match alpha_v {
        None => 0.05,
        Some(Value::Num(m)) if m.data.len() == 1 && m.data[0] > 0.0 && m.data[0] < 1.0 => m.data[0],
        _ => return rt("The 'alpha' value must be a scalar between 0 and 1."),
    };
    let freq0 = freq_in.map(|f| f.data).unwrap_or_else(|| vec![1.0; xm.data.len()]);
    let keep: Vec<bool> = xm.data.iter().zip(&freq0).map(|(a, b)| !(a.is_nan() || b.is_nan())).collect();
    let x: Vec<f64> = xm.data.iter().zip(&keep).filter(|(_, k)| **k).map(|(v, _)| *v).collect();
    let freq: Vec<f64> = freq0.iter().zip(&keep).filter(|(_, k)| **k).map(|(v, _)| *v).collect();
    let (obs, mut edges) = if let Some(Value::Num(c)) = ctrs_v {
        get_bins(it, &x, &freq, "ctrs", c)?
    } else if let Some(Value::Num(e)) = edges_v {
        get_bins(it, &x, &freq, "edges", e)?
    } else {
        let nb = match (nbins_v, exp_v) {
            (Some(Value::Num(m)), _) => m.data[0],
            (_, Some(Value::Num(e))) => e.dims.iter().copied().max().unwrap_or(0) as f64,
            _ => 10.0,
        };
        get_bins(it, &x, &freq, "nbins", &Matrix::scalar(nb))?
    };
    let mut obs = obs;
    let mut nb = obs.len();
    let mut exp: Vec<f64>;
    if let Some(ev) = exp_v {
        let bad = || rt::<()>(format!("The 'expected' value must be a vector with {nb} elements."));
        let Value::Num(e) = ev else { return bad().map(|_| None) };
        if !is_vector_dims(&e.dims) || e.data.len() != nb {
            return bad().map(|_| None);
        }
        if e.data.iter().any(|&v| v < 0.0) {
            return rt("The 'expected' values must be non-negative.");
        }
        exp = e.data.clone();
    } else {
        let sumf: f64 = freq.iter().sum();
        let (spec, mut cargs): (Value, Vec<Value>) = match cdf_v {
            None => {
                let mu = x.iter().zip(&freq).map(|(v, w)| v * w).sum::<f64>() / sumf;
                let var = x.iter().zip(&freq).map(|(v, w)| w * (v - mu) * (v - mu)).sum::<f64>() / sumf;
                if nparams.is_none() {
                    nparams = Some(2.0);
                }
                (Value::chars("normcdf"), vec![Value::num(mu), Value::num((sumf / (sumf - 1.0) * var).sqrt())])
            }
            Some(Value::Cell(c)) => {
                let rest: Vec<Value> = c.data[1..].to_vec();
                if nparams.is_none() {
                    nparams = Some(rest.len() as f64);
                }
                (c.data[0].clone(), rest)
            }
            Some(v) => (v.clone(), vec![]),
        };
        if !(spec.is_text() || matches!(spec, Value::Func(_))) {
            return rt("The 'cdf' value must be a ProbDist object, a function handle, or a cell array containing a function handle.");
        }
        let fname = crate::builtins_funfun::func_text(it, &spec);
        let interior: Vec<f64> = if edges.len() >= 2 { edges[1..edges.len() - 1].to_vec() } else { vec![] };
        let mut args = vec![row(interior)];
        args.append(&mut cargs);
        let res = match &spec {
            Value::Func(f) => it.call_value(f, args, 1),
            other => it.call_function(&other.text().unwrap_or_default(), args, 1),
        };
        let fc = match res {
            Ok(r) => r,
            Err(_) => return rt(format!("Error evaluating cdf function '{fname}'.")),
        };
        let fm = match fc.first() {
            Some(Value::Num(m)) if is_vector_dims(&m.dims) && m.data.len() == nb - 1 => m.clone(),
            _ => return rt(format!("The cdf function '{fname}' returned an incorrect number of values.")),
        };
        let total: f64 = obs.iter().sum();
        let mut prev = 0.0;
        exp = Vec::with_capacity(nb);
        for v in fm.data.iter().chain(std::iter::once(&1.0)) {
            exp.push(total * (v - prev));
            prev = *v;
        }
    }
    if exp.iter().any(|&e| e < emin) {
        let (mut i, mut j) = (0usize, exp.len() - 1);
        while i + 1 < j && (exp[i] < emin || exp[i + 1] < emin || exp[j] < emin || exp[j - 1] < emin) {
            if exp[i] < exp[j] {
                exp[i + 1] += exp[i];
                obs[i + 1] += obs[i];
                i += 1;
            } else {
                exp[j - 1] += exp[j];
                obs[j - 1] += obs[j];
                j -= 1;
            }
        }
        exp = exp[i..=j].to_vec();
        obs = obs[i..=j].to_vec();
        edges.drain(j + 1..edges.len() - 1);
        edges.drain(1..=i);
        nb = exp.len();
        if exp.iter().any(|&e| e < emin) {
            it.emit_warning("After pooling, some bins still have low expected counts.\nThe chi-square approximation may not be accurate".into(), "stats:chi2gof:LowCounts".into());
        }
    }
    let cstat: f64 = obs.iter().zip(&exp).map(|(o, e)| (o - e) * (o - e) / e).sum();
    let df = nb as f64 - 1.0 - nparams.unwrap_or(0.0);
    let (df, p) = if df > 0.0 {
        let r = it.call_function("chi2cdf", vec![Value::num(cstat), Value::num(df), Value::chars("upper")], 1)?;
        (df, r[0].scalar()?)
    } else {
        (0.0, f64::NAN)
    };
    let h = if p <= alpha { 1.0 } else { 0.0 };
    let single = is_single(xm);
    let sv = |v: f64| Value::Num(finish(Matrix::scalar(v), single));
    let srow = |v: Vec<f64>| {
        let n = v.len();
        Value::Num(finish(Matrix::from_dims(vec![1, n], v), single))
    };
    let mut out = vec![sv(h), sv(p)];
    if nargout > 2 {
        let mut st = StructArray::scalar();
        st.set("chi2stat", sv(cstat));
        st.set("df", Value::num(df));
        st.set("edges", row(edges));
        st.set("O", srow(obs));
        st.set("E", srow(exp));
        out.push(Value::Struct(st));
    }
    Ok(Some(out))
}
