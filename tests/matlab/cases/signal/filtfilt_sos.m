% filtfilt with second-order sections: each section is filtered forward and backward on its own,
% edge pad = 3 x the order of the whole cascade, steady-state start. Rounded so the last bits of sin/cos do not matter.
x = sin(0.3*(1:120)') + 0.5*cos(1.1*(1:120)'.^1.05) + 3;
[z, p, k] = butter(4, 0.3); [s, g] = zp2sos(z, p, k);
y1 = round(filtfilt(s, g, x)*1e6)/1e6;
fprintf('%.4f ', y1([1:4 60 117:120])); fprintf('\n');
[z, p, k] = ellip(6, 1, 40, 0.3); [s, g] = zp2sos(z, p, k);
y2 = round(filtfilt(s, g, x)*1e6)/1e6;
fprintf('%.4f ', y2([1:4 60 117:120])); fprintf('\n');
y3 = round(filtfilt(s, [2 3 5 7], x(1:30))*1e3)/1e3;
fprintf('%.4f ', y3([1:3 28:30])); fprintf('\n');
try, filtfilt(s, g, x(1:18)); catch e, disp(e.message), end
try, filtfilt(s, [1 2 3 4 5], x); catch e, disp(e.message), end
d = designfilt('lowpassiir', 'FilterOrder', 6, 'HalfPowerFrequency', 0.3, 'DesignMethod', 'butter');
y4 = round(filtfilt(d, x)*1e6)/1e6;
fprintf('%.4f ', y4([1:3 118:120])); fprintf('\n');
clear x z p k s g y1 y2 y3 y4 d e
