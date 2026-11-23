% thd, sfdr, snr, sinad and toi given a periodogram ('psd') or a power spectrum ('power') instead of a signal.
n = (0:1023)';
x = sin(2*pi*77.3*n/1000) + 0.1*sin(2*pi*2*77.3*n/1000+0.3) + 0.05*sin(2*pi*3*77.3*n/1000+1);
w = kaiser(1024, 38);
[P, F] = periodogram(x, w, 1024, 1000);
[S, F] = periodogram(x, w, 1024, 1000, 'power');
R = enbw(w, 1000);
[r, hp, hf] = thd(P, F, 3, 'psd');
fprintf('%.4f ', r, hp, hf); fprintf('\n');
[r, hp, hf] = thd(S, F, R, 3, 'power');
fprintf('%.4f ', r, hp, hf); fprintf('\n');
z = 0.01*cos(0.37*n.^1.3) + 0.005*sin(1.7*n) + 0.002*cos(2.9*n+0.5*n.^1.1);
[Pz, F] = periodogram(x+z+0.2, w, 1024, 1000);
[r, np] = snr(Pz, F, 3, 'psd');
fprintf('%.4f ', r, np); fprintf('\n');
[r, tp] = sinad(Pz, F, 'psd');
fprintf('%.4f ', r, tp); fprintf('\n');
[r, sp, sf] = sfdr(S, F, R, 'power');
fprintf('%.4f ', r, sp, sf); fprintf('\n');
try, sfdr(P, F, 'psd'); catch e, disp(e.message), end
try, thd(S, F, 'power'); catch e, disp(e.message), end
clear n x w P S F R r hp hf z Pz np tp sp sf e
