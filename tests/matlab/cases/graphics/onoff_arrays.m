% Arrays of matlab.lang.OnOffSwitchState values: joining, indexing, assignment, shaping
ax = axes;
a = ax.Box;
b = ax.Visible;
r = [a b]
c = [a; b]
d = [a true false]
e = [a b; b a]
s1 = r(2)
s2 = r(end)
r2 = r; r2(2) = a
r3 = r; r3(3) = b
r4 = r; r4(1) = []
r5 = r; r5(2) = 3
r6 = r'
r7 = repmat(a, 1, 3)
r8 = fliplr(r)
r9 = reshape(e, 1, 4)
q = unique([b b])
k = matlab.lang.OnOffSwitchState([true false true])
k2 = matlab.lang.OnOffSwitchState({'on', 'off'})
x = num2cell(r);
y = x{2}
t1 = double(r == 'on');
t2 = double(isequal(r, [a b]));
t3 = double(isequal(r, [false true]));
t4 = double(strcmp(r, 'on'));
t5 = double(~r);
t6 = [any(r) all(r) sum(r) numel(r) size(e, 1) size(e, 2)];
t7 = find(r);
disp(r)
disp(e)
disp(mat2str(r))
fprintf('%d %d\n', r)
fprintf('%s %s\n', r)
clear ax x
close all
