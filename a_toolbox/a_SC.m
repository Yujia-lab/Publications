function SC = a_SC(data,f);
%%% data after fft
%%% f for frequency,data for power
len = length(data);
data = reshape(data,1,len);
f = reshape(f,1,len);
SC = sum(f .* data) ./ sum(data);

