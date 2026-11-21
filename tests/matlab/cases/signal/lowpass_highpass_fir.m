% lowpass and highpass with the default steepness: the minimum-order Kaiser FIR and its zero-phase filtering.
% Only the tap counts and rounded norms are kept (the raw samples differ in the last bits of sin and cos).
x = sin(0.3*(1:400)') + 0.5*cos(1.1*(1:400)'.^1.05);
[y1, d1] = lowpass(x, 0.2);
[y2, d2] = highpass(x, 0.3);
y3 = lowpass(x, 0.2, 'Steepness', 0.9, 'StopbandAttenuation', 40);
n1 = numel(d1.Coefficients) - 1
n2 = numel(d2.Coefficients) - 1
sz = size(y3)
chk = round([norm(y1) norm(y2) norm(y3)], 3)
clear x y1 y2 y3 d1 d2
