% A tiled layout holds the axes of the tiles used so far; the figure holds the layout
f = figure;
t = tiledlayout(2, 2);
n0 = numel(f.Children);
m0 = numel(t.Children);
a1 = nexttile;
n1 = numel(f.Children);
m1 = numel(t.Children);
c1 = double(isequal(gca, a1));
a2 = nexttile;
m2 = numel(t.Children);
ord = double([t.Children(1) == a2, t.Children(2) == a1]);
a4 = nexttile(4);
m3 = numel(t.Children);
a4b = nexttile(4);
same4 = double(isequal(a4, a4b));
m4 = numel(t.Children);
g = figure;
t2 = tiledlayout(1, 1);
k0 = numel(g.Children);
l0 = numel(t2.Children);
h = figure;
t3 = tiledlayout('flow');
sz = zeros(1, 12);
for k = 1:12
    nexttile;
    sz(k) = t3.GridSize(1) * 10 + t3.GridSize(2);
end
p = figure;
t4 = tiledlayout(1, 2);
nexttile; plot(1:3);
nexttile; plot(3:-1:1);
lines = [numel(t4.Children(1).Children) numel(t4.Children(2).Children)];
q = figure;
subplot(2, 1, 1);
t5 = tiledlayout(1, 2);
cq = numel(q.Children);
clear f t a1 a2 a4 a4b g t2 h t3 p t4 q t5
close all
