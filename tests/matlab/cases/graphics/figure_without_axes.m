% A figure has no axes until something needs them; deleting the axes leaves it with none
f = figure;
n0 = numel(f.Children);
c0 = double(isempty(f.CurrentAxes));
o0 = numel(findobj(f, 'Type', 'axes'));
a = gca;
n1 = numel(f.Children);
delete(a);
n2 = numel(f.Children);
v1 = double(isvalid(a));
c2 = double(isempty(f.CurrentAxes));
g = figure; plot(1:3); clf(g);
n3 = numel(g.Children);
h = figure; ax = axes; plot(ax, 1:3); delete(ax.Children);
n4 = numel(ax.Children);
n5 = numel(h.Children);
k = figure; hold on;
n6 = numel(k.Children);
m = figure; xlim([0 1]);
n7 = numel(m.Children);
p = figure; grid on;
n8 = numel(p.Children);
q = figure; title('x');
n9 = numel(q.Children);
r = figure; colormap(r, jet);
n10 = numel(r.Children);
s = figure; cm = colormap;
n11 = numel(s.Children);
u = figure; tf = ishold;
n12 = numel(u.Children);
w = figure; axis equal;
n13 = numel(w.Children);
y = figure; ax1 = axes; ax2 = axes; delete(ax1);
n14 = numel(y.Children);
d14 = double(isequal(gca, ax2));
z = figure; ax1 = axes; ax2 = axes; delete(ax2);
d15 = double(isequal(gca, ax1));
clear f a g h ax k m p q r s u w y ax1 ax2 z cm tf
close all
