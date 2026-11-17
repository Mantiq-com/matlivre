% axes creates new axes in the figure; axes(ax) makes existing axes current
f = figure;
a1 = axes;
a2 = axes;
n1 = numel(f.Children);
s1 = double(isequal(gca, a2));
axes(a1);
s2 = double(isequal(gca, a1));
p = axes('Position', [0.1 0.2 0.5 0.4]);
pos = double(all(abs(get(p, 'Position') - [0.1 0.2 0.5 0.4]) < 1e-12));
n2 = numel(f.Children);
s3 = double(isequal(gca, p));
o1 = double(f.Children(1) == p);
o2 = double(f.Children(2) == a1);
o3 = double(f.Children(3) == a2);
g = figure;
b1 = axes(g);
b2 = axes(g);
n3 = numel(g.Children);
h = figure;
plot(1:3);
c = axes;
n4 = numel(h.Children);
n5 = numel(c.Children);
clear f a1 a2 p g b1 b2 h c
