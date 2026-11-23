% stftmag2sig with the Griffin-Lim method: zero initial phase, the stopping rule on the spectral convergence,
% the outputs t and info. Only a few samples are printed (the first iteration is the large-amplitude start).
n = (1:400)';
x = sin(0.3*n) + 0.5*sin(0.07*n.^1.2) + 0.2*cos(1.1*n);
w = hann(128, 'periodic');
S = abs(stft(x, Window=w, OverlapLength=96, FFTLength=128, FrequencyRange='centered'));
[y, t, info] = stftmag2sig(S, 128, 1000, MaxIterations=5);
fprintf('%.4f ', y(2:5), t(1:3), info.ExitFlag, info.NumIterations, info.Inconsistency); fprintf('\n');
[y, t, info] = stftmag2sig(S, 128, InconsistencyTolerance=0.5);
fprintf('%.4f ', y(2:5), info.ExitFlag, info.NumIterations, info.Inconsistency); fprintf('\n');
fprintf('%d ', size(y), size(info.ReconstructedSTFT)); fprintf('\n');
try, stftmag2sig(S, 64); catch e, disp(e.message), end
try, stftmag2sig(-S, 128); catch e, disp(e.message), end
clear n x w S y t info e
