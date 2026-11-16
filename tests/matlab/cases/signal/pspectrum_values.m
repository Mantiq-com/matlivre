% pspectrum: power spectrum and spectrogram sizes and axes (rounded so only structure is compared)
fs = 1000; n = (0:999)';
x = sin(2*pi*50*n/fs) + 0.5*sin(2*pi*120*n/fs);
[p, f] = pspectrum(x, fs);
rows = size(p, 1)
cols = size(p, 2)
df = round(f(2)*1e6)/1e6
fend = f(end)
[~, k] = max(p);
peak_hz = round(f(k))
peak_pw = round(p(k)*1e3)/1e3
[p2, f2] = pspectrum(x, fs, 'FrequencyResolution', 20);
rows2 = size(p2, 1)
step2 = f2(2)
[ps, fsg, ts] = pspectrum(x, fs, 'spectrogram');
rows_s = size(ps, 1)
cols_s = size(ps, 2)
t1 = ts(1)
dt = round((ts(2)-ts(1))*1e6)/1e6
[p3, f3] = pspectrum(x + 1i*x, fs, 'Leakage', 1);
rows3 = size(p3, 1)
f3_first = f3(1)
f3_last = f3(end)
clear x n p f p2 f2 ps fsg ts p3 f3 k
