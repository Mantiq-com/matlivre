% ifsst on a synthetic synchrosqueezed matrix: reconstruction, band and ridge inversion
nr = 9; nc = 40;
[a, b] = ndgrid(1:nr, 1:nc);
s = (sin(0.7*a + 0.3*b) + 1i*cos(0.4*a.*b/7)).*exp(-0.05*a);
f = (0:nr-1)'/(2*(nr-1));
x = ifsst(s, 16);
n1 = numel(x)
x1 = round(x(1)*1e6)/1e6
xb = ifsst(s, 16, f, [0.1 0.3]);
xb1 = round(xb(1)*1e6)/1e6
xr = ifsst(s, 16, 3*ones(nc,1));
xr1 = round(xr(1)*1e6)/1e6
clear a b s f x xb xr
