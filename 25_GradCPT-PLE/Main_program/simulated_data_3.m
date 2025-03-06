%% analyze static measurement
clear
clc


hold off;
imax = 1000;

for i = 1:imax
    i

    ple_noise(i) = -2+ (4/1000)*(i-1);

    cn = dsp.ColoredNoise('InverseFrequencyPower', ple_noise(i), 'SamplesPerFrame', 600,'NumChannels',100);
    ts = cn();

    for isub = 1:100
        x = ts(:,isub);
        SE(isub) = sampen(x,2,0.5);
        SD(isub) = std(x);
    end

    [r,p] = corrcoef(SE,SD);
    R(i) = r(2);

end

i = 751;
[AS,pos] = sort(SD);
pos_min = pos(1:50);
pos_max = pos(51:end);

ts_min = ts(:,pos_min);
ts_max = ts(:,pos_max);

psd_min = pwelch(ts_min,[],[],[],fs);
psd_max = pwelch(ts_max,[],[],[],fs);

mean_psd_min = mean(psd_min,2);
mean_psd_max = mean(psd_max,2);
mean_sd_min = mean(SD(pos_min));
mean_sd_max = mean(SD(pos_max));
mean_se_min = mean(SE(pos_min));
mean_se_max = mean(SE(pos_max));

for i = 1:129
    x1 = psd_min(i,:);
    x2 = psd_max(i,:);
    [h,p,ci,stats] = ttest2(x2,x1);
    P_psd(i) = p;
    T_psd(i) = stats.tstat;
end

[psd1, freq] = pwelch(ts(:,9), [], [],[], fs);
[psd2, freq] = pwelch(ts(:,48), [], [], [], fs);


[psd, freq] = periodogram(ts(:,2), [], [], 1);

coeffs = polyfit(log10(freq), log10(psd), 1);
PLE = -coeffs(1);
y = coeffs(2) + coeffs(1)*log10(freq2);
