% Multiple assignment into comma-separated lists: output k goes to element k, in
% column-major order. Results are printed as text (the command-window display of cells and
% struct arrays is covered elsewhere).
[c{1:2}] = deal(1, 2);
fprintf('c: %s\n', mat2str([c{:}]));
[sz{1:2}] = size(zeros(2, 3));
fprintf('sz: %s\n', mat2str([sz{:}]));
d = cell(1, 3); [d{:}] = deal('a', 'b', 'c');
fprintf('d: %s\n', [d{:}]);
[a, g{2:3}, b] = deal(1, 2, 3, 4);
fprintf('a=%d b=%d g: %s, isempty(g{1})=%d\n', a, b, mat2str([g{2:3}]), isempty(g{1}));
h = cell(2, 2); [h{:, 2}] = deal(8, 9);
fprintf('h: %s %s\n', mat2str(cellfun(@numel, h)), mat2str([h{:, 2}]));
[k{2, 1:2}] = deal(3, 4);
fprintf('k: %s %s\n', mat2str(size(k)), mat2str([k{2, :}]));
m = {1}; [m{1:3}] = deal('x', 'y', 'z');
fprintf('m: %s\n', [m{:}]);
p = {1, 2, 3}; [p{[3 1]}] = deal('p', 'q');
fprintf('p: %s %d %s\n', p{1}, p{2}, p{3});
q = {1, 2, 3}; [q{logical([1 0 1])}] = deal(5, 6);
fprintf('q: %s\n', mat2str([q{:}]));
r = {1, 2}; [~, r{:}] = deal(7, 8, 9);
fprintf('r: %s\n', mat2str([r{:}]));
t.c = cell(1, 2); [t.c{:}] = deal(5, 6);
fprintf('t.c: %s\n', mat2str([t.c{:}]));
z = {}; [z{:}] = deal();
fprintf('z: %s\n', mat2str(size(z)));

% Struct arrays: [s.f] and [s(idx).f].
s = struct('f', {1, 2, 3});
[s.f] = deal(10, 20, 30);
fprintf('s.f: %s\n', mat2str([s.f]));
[s.f] = deal(7);
fprintf('s.f: %s\n', mat2str([s.f]));
[s(2:3).f] = deal('u', 'v');
fprintf('s.f: %d %s %s\n', s.f);
u = struct('f', {1, 2}); [u(3:4).f] = deal(5, 6);
fprintf('u: %s %s\n', mat2str(size(u)), mat2str([u.f]));
v = struct('f', {1, 2}); [v.f] = size(ones(4, 5));
fprintf('v.f: %s\n', mat2str([v.f]));
e = struct('f', {}); [e.f] = deal(1);
fprintf('e: %s\n', mat2str(size(e)));
[w.f] = deal(11);
fprintf('w: %s %d\n', mat2str(size(w)), w.f);

% Errors.
try, [n1{1:2}] = numel(1); catch ME, fprintf('%s | %s\n', ME.identifier, ME.message); end
try, n2 = {}; [n2{1:2}] = 5; catch ME, fprintf('%s | %s\n', ME.identifier, ME.message); end
try, n3 = struct('f', {1, 2}); [n3.f] = 5; catch ME, fprintf('%s | %s\n', ME.identifier, ME.message); end
try, n4 = cell(1, 2); [n4{:}] = deal(1, 2, 3); catch ME, fprintf('%s | %s\n', ME.identifier, ME.message); end
try, n5 = 5; [n5{1:2}] = deal(1, 2); catch ME, fprintf('%s | %s\n', ME.identifier, ME.message); end
try, n6 = 5; n6{1} = 2; catch ME, fprintf('%s | %s\n', ME.identifier, ME.message); end
try, n7 = {}; n7{1:2} = 5; catch ME, fprintf('%s | %s\n', ME.identifier, ME.message); end
try, n8 = {1, 2}; n8{:} = 5; catch ME, fprintf('%s | %s\n', ME.identifier, ME.message); end
try, n9 = struct('f', {1, 2}); n9.f = 5; catch ME, fprintf('%s | %s\n', ME.identifier, ME.message); end
try, n10 = struct('f', {1, 2}); [n10.f.g] = deal(1, 2); catch ME, fprintf('%s | %s\n', ME.identifier, ME.message); end
try, [n11{1:2}] = undefinedFunctionForCsList(); catch ME, fprintf('%s | %s\n', ME.identifier, ME.message); end
try, n12{undefinedIndexForCsList} = undefinedFunctionForCsList(); catch ME, fprintf('%s | %s\n', ME.identifier, ME.message); end
try, n14 = {}; n14{:} = 5; catch ME, fprintf('%s | %s\n', ME.identifier, ME.message); end
try, n15 = 5; n15{1:2} = 3; catch ME, fprintf('%s | %s\n', ME.identifier, ME.message); end
try, [n16, n17] = 5; catch ME, fprintf('%s | %s\n', ME.identifier, ME.message); end
n13 = {1}; n13{:} = 5;
fprintf('n13: %d\n', n13{1});
n5, n6, n15

clear ME c sz d g h k m p q r t z s u v e w n1 n2 n3 n4 n7 n8 n9 n10 n11 n12 n13 n14
