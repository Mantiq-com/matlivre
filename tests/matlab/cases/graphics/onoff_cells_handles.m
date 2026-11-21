% OnOffSwitchState values inside cells, upper/lower, isequal on handles, anonymous functions
ax = axes;
v = ax.Visible;
r1 = double(any(strcmp({v}, 'on')));
r2 = double(strcmp({v, 'x'}, 'on'));
r3 = double(isequal({v}, {'on'}));
r4 = double(isequal(gca, gca));
a1 = gca;
c = get(gcf, 'Children');
r5 = double(isequal(c(1), a1));
r6 = double(isequal(v, upper(v)));
disp(class(upper(v)))
disp(class(lower(v)))
f = @(x) fprintf('%d\n', x);
f(3)
g = @(x) x + 1;
g(2)
clear ax v a1 c f g
close all
