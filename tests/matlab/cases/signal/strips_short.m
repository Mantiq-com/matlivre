% strips with a single strip (signal shorter than the strip) and a constant signal.
t = (0:29)';
strips(round(sin(0.4*t)*1e6)/1e6, 100);
figure;
strips(ones(400,1), 100);
figure;
strips(round(sin(0.3*(1:25)')*1e6)/1e6, 10);
clear t
