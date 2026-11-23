% yulewalk: recursive filter design by a least-squares fit to a magnitude response
f = [0 .55 .6 .65 1];
m = [1 1 .5 0 0];
[b, a] = yulewalk(8, f, m);
sz = [size(b) size(a)]
a1 = a(1)
bq = round(b*1e6)
aq = round(a*1e6)
stable = double(all(abs(roots(a)) < 1))
h = freqz(b, a, [0 0.3 0.9]*pi);
hq = round(abs(h)*1e4)
[b2, a2] = yulewalk(6, [0 .2 .2 .5 .5 1]', [0 0 1 1 0 0]', 256, 8);
sz2 = [size(b2) size(a2)]
b2q = round(b2*1e6)
a2q = round(a2*1e6)
b3 = yulewalk(3, single([0 1]), single([1 0.25]));
is_double3 = double(isa(b3, 'double'))
b3q = round(b3*1e6)
try
    yulewalk(4, [0 .6 .4 1], [1 1 0 0]);
catch err
    msg1_ok = double(strcmp(err.message, 'Frequencies must be non-decreasing.'))
end
try
    yulewalk(6, [0 0 .5 1], [0 1 1 1]);
catch err
    msg2_ok = double(strcmp(err.message, 'Too abrupt an amplitude change near end of frequency interval.'))
end
clear f m b a h b2 a2 b3 err
