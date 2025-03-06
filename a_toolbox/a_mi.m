function MI = a_mi(x)

% calculate the modulation index

mean_x = x ./ sum(x);

hp = -sum(mean_x.*log10(mean_x));
MI = (log10(length(mean_x))-hp)/log10(length(mean_x));