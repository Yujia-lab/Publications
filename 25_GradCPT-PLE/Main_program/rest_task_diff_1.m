%% analyze static measurement
clear
clc

rest_path = 'E:\GradCPT\ROISignals_filt_regressed\ses-1';
task_path = 'E:\GradCPT\ROISignals_filt_regressed\ses-2';

rest_file = dir(rest_path);
rest_file(1:2) = [];
task_file = dir(task_path);
task_file(1:2) = [];

rest_path2 = 'E:\GradCPT\network_signals_regressed\ses-1';
task_path2 = 'E:\GradCPT\network_signals_regressed\ses-2';

rest_file2 = dir(rest_path2);
rest_file2(1:2) = [];
task_file2 = dir(task_path2);
task_file2(1:2) = [];

TR=0.8;%TR
nfft=512;
window=hamming(32);
overlap=16;
%  len=(nfft/2+1)*2;
fs=1/TR;


for isubj = 1:length(rest_file)
    isubj
    load(fullfile(rest_path,rest_file(isubj).name));
    rest_signal = roisignals_filt_rest;
    load(fullfile(task_path,task_file(isubj).name));
    task_signal = roisignals_filt_task;

     load(fullfile(rest_path2,rest_file2(isubj).name));
     rest_signal2 = roisignals(:,4:end)';
     load(fullfile(task_path2,task_file2(isubj).name));
     task_signal2 = roisignals(:,4:end)';

    for iroi = 1:size(rest_signal,2)
        x1 = rest_signal(:,iroi);
        x2 = task_signal(:,iroi);

         x11 = rest_signal2(:,iroi);
         x22 = task_signal2(:,iroi);

        sd_rest(iroi,isubj) = std(x1);
        sd_task(iroi,isubj) = std(x2);

        [ple_rest(iroi,isubj),psd_rest2(:,iroi,isubj),freq] = ple(x11,fs,[0.01,0.1],0);
        [ple_task(iroi,isubj),psd_task2(:,iroi,isubj)] = ple(x22,fs,[0.01,0.1],0);

        sampen_rest(iroi,isubj) = sampen(x1,2,0.5);
        sampen_task(iroi,isubj) = sampen(x2,2,0.5);

        [psd_rest(:,iroi,isubj),f] = pwelch(x11,[],[],nfft,fs);
        [psd_task(:,iroi,isubj),f] = pwelch(x22,[],[],nfft,fs);
    end
end



for i = 1:7
    sd_rest1 = sd_rest(i,:);
    sampen_rest1 = sampen_rest(i,:);
    [r,p] = corrcoef(sd_rest1,sampen_rest1);
    R(i) = r(2);
end

mean_sd_rest = mean(sd_rest,2);
mean_sd_task = mean(sd_task,2);

mean_sampen_rest = mean(sampen_rest,2);
mean_sampen_task = mean(sampen_task,2);

mean_ple_rest = mean(ple_rest,2);
mean_ple_task = mean(ple_task,2);

path = 'E:\GradCPT\visualization';
yeo_visualize(mean_sd_rest(2:8),path,'mean_sd_rest')
yeo_visualize(mean_sd_task(2:8),path,'mean_sd_task')
yeo_visualize(mean_sampen_rest(2:8),path,'mean_sampen_rest')
yeo_visualize(mean_sampen_task(2:8),path,'mean_sampen_task')
yeo_visualize(mean_ple_rest(2:8),path,'mean_ple_rest')
yeo_visualize(mean_ple_task(2:8),path,'mean_ple_task')


for iroi = 1:8
    x1 = sd_rest(iroi,:);
    x2 = sd_task(iroi,:);
    [h,p,ci,stats] = ttest(x1,x2);
    P_sd(iroi) = p;
    T_sd(iroi) = stats.tstat;

    x1 = sampen_rest(iroi,:);
    x2 = sampen_task(iroi,:);
    [h,p,ci,stats] = ttest(x1,x2);
    P_sampen(iroi) = p;
    T_sampen(iroi) = stats.tstat;

    x1 = ple_rest(iroi,:);
    x2 = ple_task(iroi,:);
    [h,p,ci,stats] = ttest(x1,x2);
    P_ple(iroi) = p;
    T_ple(iroi) = stats.tstat;
end

T_sd_FDR = a_multicorrect(T_sd,P_sd,'FDR');
T_sampen_FDR = a_multicorrect(T_sampen,P_sampen,'FDR');
T_ple_FDR = a_multicorrect(T_ple,P_ple,'FDR');









