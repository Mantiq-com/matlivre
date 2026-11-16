% orderspectrum: average order spectrum of a tachometer-referenced signal
fs = 1000; n = 1024; t = (0:n-1)'/fs;
rpm = 600 + 300*t + 40*sin(2*pi*1.5*t);
ph = 2*pi*cumtrapz(t, rpm/60);
x = cos(ph) + 0.5*cos(2*ph) + 0.1*sin(2*pi*90*t);
[spec, order] = orderspectrum(x, fs, rpm);
npts = numel(spec)
step = round(order(2)*1e4)/1e4
[~, k] = max(spec);
peak_order = round(order(k))
peak_amp = round(spec(k)*100)/100
clear t rpm ph x spec order k
