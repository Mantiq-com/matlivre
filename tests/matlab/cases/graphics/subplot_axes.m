% subplot makes or reuses one axes per call and deletes the axes it overlaps
near = @(p, q) double(all(abs(p - q) < 1e-6));
f = figure;
a = subplot(2,1,1);
n1 = numel(f.Children);
b = subplot(2,1,2);
n2 = numel(f.Children);
a2 = subplot(2,1,1);
n3 = numel(f.Children);
same = double(isequal(a, a2));
cur = double(isequal(gca, a));
front = double(isequal(f.Children(1), b));
p1 = near(a.Position, [0.13 0.583837 0.775 0.341163]);
p2 = near(b.Position, [0.13 0.11 0.775 0.341163]);
g = figure;
a = subplot(2,1,1);
b = subplot(2,1,2);
c = subplot(1,2,1);
ov1 = numel(g.Children);
ov2 = [double(isvalid(a)) double(isvalid(b)) double(isvalid(c))];
pc = near(c.Position, [0.13 0.11 0.334659 0.815]);
h = figure;
a = subplot(2,2,1); b = subplot(2,2,2); c = subplot(2,2,3); d = subplot(2,2,4);
e = subplot(1,2,2);
ov3 = numel(h.Children);
ov4 = [double(isvalid(a)) double(isvalid(b)) double(isvalid(c)) double(isvalid(d))];
pe = near(e.Position, [0.570341 0.11 0.334659 0.815]);
k = figure;
s1 = subplot(2,2,[1 2]);
s2 = subplot(2,2,3);
s3 = subplot(2,2,4);
sp1 = near(s1.Position, [0.13 0.583837 0.775 0.341163]);
sp2 = near(s2.Position, [0.13 0.11 0.334659 0.341163]);
sp3 = near(s3.Position, [0.570341 0.11 0.334659 0.341163]);
sk = numel(k.Children);
m = figure;
t1 = subplot(3,3,[2 9]);
st = near(t1.Position, [0.410797 0.11 0.494203 0.815]);
q = figure;
w1 = subplot(221);
w2 = subplot('Position', [0.2 0.3 0.4 0.5]);
wp = near(w2.Position, [0.2 0.3 0.4 0.5]);
wn = numel(q.Children);
r = figure;
x1 = subplot(2,1,1);
x2 = subplot(2,1,1, 'replace');
xv = [double(isvalid(x1)) double(isvalid(x2))];
xn = numel(r.Children);
u = figure;
y1 = subplot(2,1,1);
y2 = subplot(2,1,2);
subplot(y1);
yc = double(isequal(gca, y1));
yn = numel(u.Children);
v = figure;
subplot(2,1,1); plot(1:3);
subplot(2,1,2); plot(3:-1:1);
vl = [numel(v.Children(1).Children) numel(v.Children(2).Children)];
z = figure;
plot(1:3);
z1 = gca;
z2 = subplot(1,1,1);
zs = double(isequal(z1, z2));
zl = numel(z2.Children);
zn = numel(z.Children);
y = figure;
plot(1:3);
y3 = gca;
subplot(2,1,1);
yd = double(isvalid(y3));
yk = numel(y.Children);
try
    subplot(2,2,5)
catch err
    disp(err.message)
end
try
    subplot(2,2,0)
catch err
    disp(err.message)
end
clear near f a b a2 g c h d e k s1 s2 s3 m t1 q w1 w2 r x1 x2 u y1 y2 v z z1 z2 y y3 err
close all
