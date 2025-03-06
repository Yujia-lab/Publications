function [lzc] = LZC(x)
% LZC calculates the Lempel-Ziv complexity of a non-binary signal by first binarizing it according
% to its median, then applying the LZ algorithm, and finally normalizing by its length.
% Adapted from the wikipedia page by Yasir Çatal aka Duodenum
% Set parameters
n = length(x);
b = n/log2(n);
threshold = median(x);

% Binary series of 'ones' (higher values than median) and 'zeros' (lower values than median) 
bin = x; % Initialize binarized vector
bin(x >= threshold) = 1;
bin(x < threshold) = 0;

% for j = 1:n
%     if x(j) >= threshold
%    	    x(j) = 1;
% 	else
%    	    x(j) = 0;
%     end
% end

%% LZC algorithm, check the wiki page of LZC for the meaning of variables
i = 0;
c = 1;
u = 1;
v = 1;
vmax = v;

while u + v <= n
    if bin(i + v) == bin(u + v)
        v = v + 1;
    else
        vmax = max(v, vmax);
        i = i + 1;
        if i == u
            c = c + 1;
            u = u + vmax;
            v = 1;
            i = 0;
            vmax = v;
        else
            v = 1;
        end
    end
end
% if v ~= 1
%     c = c + 1;
% end
lzc  = c / b;
end
            
            
            
            
            
            


