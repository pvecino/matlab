function [x, n] = expon(b, n0, L)
n = n0:(n0 + L -1);
x = b.^n;
end