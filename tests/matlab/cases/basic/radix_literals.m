% Hex and binary literals are integers: the smallest unsigned class that holds them, or the one a suffix names
a = 0x1F;
b = 0xFF;
c = 0x100;
d = 0xFFFF;
e = 0x10000;
f = 0xFFFFFFFF;
g = 0x100000000;
h = 0b101;
i1 = 0xFFu8;
i2 = 0x7Fs16;
i3 = 0x1Fu16;
i4 = 0b1u8;
i5 = 0x7Fs8;
i6 = 0xFFFFs32;
i7 = 0b1111s64;
i8 = 0xABCdef;
i9 = 0XA;
fn = @(x) 0x1F + x;
r1 = fn(2);
r2 = fn(250);
r3 = fn(1.5);
t1 = func2str(fn);
t2 = func2str(@(y) 0b101*y + 0xFFu8);
v = [0x10 0b11 5];
w = 0x10 + 0x01;
cls = {class(a), class(b), class(c), class(e), class(g), class(h), class(i2), class(i5), class(v), class(w)};
fprintf('%s\n', strjoin(cls, ' '));
q = 1:0x4;
z = zeros(0x2, 0b11);
sz = size(z);
% a 0 followed by an identifier is still a product-free syntax error elsewhere; plain digits are unchanged
p = 10;
p2 = p*0 + 0.5e1 + 1e0;
clear cls
