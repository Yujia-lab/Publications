function x_angle = a_Phase(data)

%%% calculate phase for timeseries

x = hilbert(data); % hilbert trans
x_angle = angle(x);
% delete first and last 10 points
% x_angle((end-9):end) = [];
% x_angle(1:10) = [];
x_angle = x_angle * 180 / pi; 
