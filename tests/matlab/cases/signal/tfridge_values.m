% tfridge: ridge extraction by penalized dynamic programming
nr = 40; nc = 60;
[a, b] = ndgrid(1:nr, 1:nc);
s = exp(-((a-10-0.3*b-4*sin(b/6)).^2)/18) + 0.7*exp(-((a-30+0.1*b).^2)/8) + 0.1*abs(sin(a.*b*0.37));
f = (0:nr-1)'*5;
[fr, ir, lr] = tfridge(s, f);
sz = size(fr)
first = fr(1)
last = ir(end)
lin_ok = double(all(lr == sub2ind(size(s), ir, (1:nc)')))
[fr2, ir2] = tfridge(s, f, 5, 'NumRidges', 2);
sz2 = size(ir2)
r2 = ir2(1,2)
clear a b s f fr ir lr fr2 ir2
