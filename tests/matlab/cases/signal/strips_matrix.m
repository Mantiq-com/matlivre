% strips of a matrix: the columns share the strips, one midrange offset for all, and SCALE.
t = (0:949)';
x = round((sin(2*pi*t/100) + 0.3*cos(2*pi*t/37))*1e6)/1e6;
strips([x 2*x+3 x.^2], 2, 100, 0.5);
clear t x
