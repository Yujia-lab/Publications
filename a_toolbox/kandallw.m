function w = kandallw(x)

[n,k] = size(x);
SR = sum(x,2); SRBAR = mean(SR);
S = sum(SR.^2) - n*SRBAR^2;
w = 12*S/k^2/(n^3-n);