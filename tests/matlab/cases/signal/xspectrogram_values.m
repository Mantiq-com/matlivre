% xspectrogram: cross-spectrogram sizes, axes and magnitudes
fs = 1000; n = (0:799)';
x = sin(2*pi*50*n/fs) + 0.3*sin(2*pi*210*n/fs);
y = filter([0.5 0.4], [1 -0.3], x);
[s, f, t, c] = xspectrogram(x, y, 64, 32, 128, fs);
rows = size(s, 1)
cols = size(s, 2)
df = f(2)
t1 = t(1)
peak_row = find(s(:,5) == max(s(:,5)))
mag_matches = double(max(max(abs(abs(c) - s))) < 1e-12)
[s2, f2, t2] = xspectrogram(x, y);
rows2 = size(s2, 1)
cols2 = size(s2, 2)
clear x y n s f t c s2 f2 t2
