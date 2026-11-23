% strips with a strip length in samples; the last strip is padded with NaN; labels are start samples.
t = (0:949)';
x = round((sin(2*pi*t/100) + 0.3*cos(2*pi*t/37))*1e6)/1e6;
strips(x, 200);
clear t x
