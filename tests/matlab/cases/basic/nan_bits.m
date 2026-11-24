% The sign bit of a NaN is not part of what is compared: NaNs made in different ways are the same value
a = NaN;
b = 0/0;
c = [1 NaN 3];
d = single(NaN);
e = Inf - Inf;
f = nan(2, 2);
h = -NaN;
k = [NaN Inf -Inf];
n = NaN * 0;
