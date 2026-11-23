% thd, sfdr and toi on noise-free tones (all harmonics present). Printed with fixed digits so the
% last bits of sin and cos do not matter.
n = (0:1023)';
x = sin(2*pi*77.3*n/1000) + 0.1*sin(2*pi*2*77.3*n/1000+0.3) + 0.05*sin(2*pi*3*77.3*n/1000+1);
[r, hp, hf] = thd(x, 1000, 3);
fprintf('%.4f ', r, hp, hf); fprintf('\n');
[r, sp, sf] = sfdr(x, 1000);
fprintf('%.4f ', r, sp, sf); fprintf('\n');
[r, sp, sf] = sfdr(x, 1000, 80);
fprintf('%.4f ', r, sp, sf); fprintf('\n');
t = (0:2047)';
y = 0.01*sin(2*pi*200*t/2000) + 0.02*sin(2*pi*250*t/2000);
y = y + 1e-2*y.^3;
[p, fp, ff, ip, fi] = toi(y, 2000);
fprintf('%.4f ', p, fp, ff, ip, fi); fprintf('\n');
try, thd(x, 1000, 1); catch e, disp(e.message), end
try, sfdr(x, 1000, -1); catch e, disp(e.message), end
try, toi([x x], 1000); catch e, disp(e.message), end
z = 0.01*cos(0.37*n.^1.3) + 0.005*sin(1.7*n) + 0.002*cos(2.9*n+0.5*n.^1.1);
w = x + z + 0.2;
[r, np] = snr(w, 1000, 3);
fprintf('%.4f ', r, np); fprintf('\n');
[r, tp] = sinad(w, 1000);
fprintf('%.4f ', r, tp); fprintf('\n');
r = snr(x, z);
fprintf('%.4f\n', r);
clear n x z w np tp r hp hf sp sf t y p fp ff ip fi e
