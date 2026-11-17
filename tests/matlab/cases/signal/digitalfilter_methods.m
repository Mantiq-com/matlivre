% digitalFilter object functions
x = sin(0.3*(1:40)') + 0.5*cos(1.1*(1:40)'.^1.05);
d1 = designfilt('lowpassfir','FilterOrder',20,'PassbandFrequency',0.2,'StopbandFrequency',0.3);
d2 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',0.3);
isfir1 = isfir(d1)
isfir2 = isfir(d2)
ord1 = filtord(d1)
stab = double(isstable(d2))
y1 = round(sum(filter(d1, x))*1e6)/1e6
y2 = round(sum(filter(d2, x))*1e6)/1e6
[b, a] = tf(d2);
nb = numel(b)
a2 = round(a(2)*1e6)/1e6
[h, w] = freqz(d1, 64);
nh = numel(h)
w2 = round(w(2)*1e6)/1e6
same = double(d1 == d1)
dif = double(d1 ~= d2)
clear x d1 d2 b a h w y1 y2
