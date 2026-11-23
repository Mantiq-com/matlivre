% fillgaps: forward pass (left to right) and backward pass (right to left) AR predictions, mean removed,
% weighted k/(L+1) across the gap; fixed order and automatic order (AIC), with a search length.
n = (1:60)';
x = sin(0.31*n) + 0.5*cos(0.77*n.^1.1) + 0.3*sin(1.9*n+0.4) + 0.001*n;
y = x; y(20:22) = NaN; y(40:41) = NaN;
a = fillgaps(y, [], 3);
fprintf('%.6f ', a([20:22 40:41])); fprintf('\n');
z = x + 0.05*sin(12345.678*n.^2);
z(30:32) = NaN;
b = fillgaps(z);
fprintf('%.6f ', b(30:32)); fprintf('\n');
c = fillgaps(z, 20);
fprintf('%.6f ', c(30:32)); fprintf('\n');
d = fillgaps(z.', [], 4);
fprintf('%.6f ', d(30:32)); fprintf('\n');
try, fillgaps(z, 8, 20); catch e, disp(e.message), end
try, fillgaps(z, [], 0); catch e, disp(e.message), end
clear n x y z a b c d e
