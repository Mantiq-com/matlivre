% modalfit (least-squares complex exponential) on a two-mode receptance FRF
fs = 1000; N = 512; f = (0:N/2)'*fs/N; w = 2*pi*f;
H = zeros(size(f));
H = H + 1./((2*pi*50)^2 - w.^2 + 2i*0.02*(2*pi*50)*w);
H = H + 0.5./((2*pi*120)^2 - w.^2 + 2i*0.03*(2*pi*120)*w);
[fn, dr, ms, ofrf] = modalfit(H, f, fs, 3);
fn1 = round(fn(1)*1e3)/1e3
fn2 = round(fn(2)*1e3)/1e3
fn3_is_nan = double(isnan(fn(3)))
dr1 = round(dr(1)*1e5)/1e5
dr2 = round(dr(2)*1e5)/1e5
nms = numel(ms)
nof = numel(ofrf)
clear f w H fn dr ms ofrf
