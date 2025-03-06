function tsout = a_FIRfilter(ts, band, Fs)

%%ts填你要滤的一列或者一行数据
%%band填要滤的频段，比如说[0.02,0.05]
%%FS填扫描频率，一般是0.5

DEV = [0.05 0.01 0.05];
A = [0 1 0];
F = [band(1)-0.01, band(1), band(2), band(2)+0.01];
[N, Fo, Ao, W] = firpmord(F, A, DEV, Fs);
if(mod(N, 2) == 1) 
    N = N + 1; % it has to be even -> b is odd -> delay is integer = N/2
end
psess.filter.b = firpm(N, Fo, Ao, W);
psess.filter.N = N;
         
tsnorm = ts - mean(ts);             % removing DC for optimal bandpass filtering
T = length(tsnorm);
tsout = conv(psess.filter.b, tsnorm);  % produces a T+N-1 length signal
tsout = flipud(conv(psess.filter.b, flipud(tsout)));
if(0) % delay compensation, normal filter [not implemented]
    tsout = tsout((psess.filter.N/2 + 1):(end - psess.filter.N/2));
else
    tsout = tsout(length(psess.filter.b):end);
    tsout(T+1:end) = [];
end
if(0) % option for filter transient removal
    tsout = tsout((psess.filter.N + 1):end); % removing N filter transient and N/2 filter tail of factual data
end
