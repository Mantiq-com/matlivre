% strips: a long signal plotted in horizontal strips (one line with NaN separators), the mean removed,
% the strip spacing 1.25 times the range of the data, tick labels at the strip start times.
t = (0:999)';
x = round((sin(2*pi*t/100) + 0.3*cos(2*pi*t/37))*1e6)/1e6;
strips(x, 2, 100);
clear t x
