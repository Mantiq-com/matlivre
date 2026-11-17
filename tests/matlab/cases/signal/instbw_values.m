% instfreq (spectrogram moment method) and instbw on a linear chirp
fs = 1000; n = (0:599)';
x = cos(2*pi*(40*n/fs + 30*(n/fs).^2));
[f1, t1] = instfreq(x, fs);
npts = numel(f1)
tfirst = t1(1)
mid = round(f1(10))
[bw, tb] = instbw(x, fs);
nbw = numel(bw)
bwmid = round(bw(10))
clear x n f1 t1 bw tb
