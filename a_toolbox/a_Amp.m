function x_amp = a_Amp(xn)
%%% Coding by Aoyujia 2020.7.20
%%% See manu for more instruction °´ÁÐ¼ÆËãÏàÎ»

xn1 = hilbert(xn);
xr = real(xn1);
xi = imag(xn1);
x_amp = sqrt(xr.^2+xi.^2);
x_amp((end-9):end) = [];
x_amp(1:10) = [];