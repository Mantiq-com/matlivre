% savefig and openfig keep where the axes are
near = @(p, q) double(all(abs(p - q) < 1e-6));
fn = fullfile(tempdir, 'livre_case_axes_positions.fig');
f = figure;
subplot(2, 1, 1); plot(1:3);
subplot(2, 1, 2); plot(3:-1:1);
savefig(f, fn);
close(f);
g = openfig(fn);
n1 = numel(g.Children);
p1 = near(g.Children(2).Position, [0.13 0.583837 0.775 0.341163]);
p2 = near(g.Children(1).Position, [0.13 0.11 0.775 0.341163]);
d1 = [g.Children(2).Children.YData; g.Children(1).Children.YData];
close(g);
f = figure;
a = axes('Position', [0.1 0.2 0.5 0.4]);
plot(a, 1:3);
savefig(f, fn);
close(f);
g = openfig(fn);
n2 = numel(g.Children);
p3 = near(g.Children(1).Position, [0.1 0.2 0.5 0.4]);
close(g);
f = figure;
a = axes;
b = axes('Position', [0.2 0.2 0.3 0.3]);
plot(a, 1:3);
plot(b, 3:-1:1);
savefig(f, fn);
close(f);
g = openfig(fn);
n3 = numel(g.Children);
p4 = near(g.Children(2).Position, [0.13 0.11 0.775 0.815]);
p5 = near(g.Children(1).Position, [0.2 0.2 0.3 0.3]);
close(g);
delete(fn);
clear near fn f a b g
close all
