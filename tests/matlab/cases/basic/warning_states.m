% warning states: the state a change returns, restoring it, 'all', and identifiers
A = [1 1; 1 1]; b = [1; 2];
show = @(w) disp(strjoin(arrayfun(@(e) [e.identifier '=' e.state], w, 'UniformOutput', false), ' '));
s = warning('off', 'MATLAB:singularMatrix');
show(s)
show(warning('query', 'MATLAB:singularMatrix'))
warning('query', 'MATLAB:singularMatrix')
warning(s);
show(warning('query', 'MATLAB:singularMatrix'))
a = warning('off', 'all');
n_a = numel(a);
show(warning('query', 'my:id'))
x = A \ b;
w1 = warning('on', 'all');
n_w1 = numel(w1);
warning(a);
% a disabled warning is not shown but is remembered by lastwarn
warning('off', 'my:id');
lastwarn('');
warning('my:id', 'hello %d', 3);
[m1, id1] = lastwarn;
disp([m1 ' / ' id1])
% 'error' turns the warning into an error
warning('error', 'my:id');
try
    warning('my:id', 'stop %s', 'here');
catch err
    disp([err.identifier ' / ' err.message])
end
% the saved states, in the order they were set
warning('on', 'all');
warning('off', 'a:b');
warning('error', 'c:d');
w = warning('query', 'all');
n = numel(w);
show(w)
w2 = warning('off', 'all');
n2 = numel(w2);
w3 = warning('on', 'all');
n3 = numel(w3);
warning(w2);
w4 = warning('query', 'all');
same = double(isequal(w4, w));
% built-in warnings have identifiers (a disabled one is still remembered by lastwarn)
warning('on', 'all');
warning('off', 'MATLAB:singularMatrix');
lastwarn('');
y = inv(A);
[m2, id2] = lastwarn;
disp([m2 ' / ' id2])
% state of the backtrace mode
show(warning('off', 'backtrace'))
warning('on', 'backtrace');
warning('on', 'all');
clear show A b s a w1 m1 id1 err w w2 w3 w4 m2 id2
