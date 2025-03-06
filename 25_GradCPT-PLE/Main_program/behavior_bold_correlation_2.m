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

% psd shape
psd_all1 = mean(psd_rest2,3);
psd_all2 = mean(psd_all1,2);

path = 'E:\GradCPT\Behavior';
subj_file = dir(path);
subj_file(1:2) = [];

% 插值
for isubj = 1:length(subj_file)
    load(fullfile(path,subj_file(isubj).name));
    result2 = result.final;
    data = squeeze(struct2cell(result2));

    data(2,:) = [];

    RT = cell2mat(data(7,:));
    cellArray = data(7,:); 
    vector = NaN(1, 400);
    
    for i = 1:length(cellArray)
        if ~isempty(cellArray{i})
            vector(i) = cellArray{i};
        end
    end

    nonNanIndices = find(~isnan(vector));
    nonNanValues = vector(nonNanIndices);

    med = median(nonNanValues);
    vector2 = vector;
    vector2(isnan(vector2)) = med;
    
    % 对整个向量范围进行插值
    interpVector = interp1(nonNanIndices, nonNanValues, 1:length(vector), 'spline', 'extrap');
    RT_series(:,isubj) = interpVector;
    %RT_series(:,isubj) = vector2;

    Acc(:,isubj) = cell2mat(data(8,:))';
end

RT_series(1:2,:) = [];
RT_series(end-1:end,:) = [];

Acc(1:2,:) = [];
Acc(end-1:end,:) = [];

RT_sd = std(RT_series,0,1)';

for i = 1:49
    RT_lzc(i) = LZC(RT_series(:,i));
    [rt_ple(i),rt_psd(:,i)] = ple(RT_series(:,i),1/1.2,[0.01,0.1],0);
    ACC_com(i) = length(find(Acc(:,i)==3))/100;
    ACC_om(i) = length(find(Acc(:,i)==4))/100;
end

[R1,P1] = corr(sampen_task',RT_sd,'type','Pearson');
[R2,P2] = corr(sampen_rest',RT_sd,'type','Pearson');

R_sampen = R1;
P_sampen = P1;

R_FDR1 = a_multicorrect(R1(2:8),P1(2:8),'FDR');
R_FDR2 = a_multicorrect(R2(2:8),P2(2:8),'FDR');

path = 'E:\GradCPT\visualization';
yeo_visualize(R_FDR1,path,'correlation_sampen_task')
yeo_visualize(R_FDR2,path,'correlation_sampen_rest')

[R1,P1] = corr(sd_task',RT_sd,'type','Pearson');
[R2,P2] = corr(sd_rest',RT_sd,'type','Pearson');

R_sd = R1;
P_sd = P1;

R_FDR1 = a_multicorrect(R1(2:8),P1(2:8),'FDR');
R_FDR2 = a_multicorrect(R2(2:8),P2(2:8),'FDR');

path = 'E:\GradCPT\visualization';
yeo_visualize(R_FDR1,path,'correlation_sd_task')
yeo_visualize(R_FDR2,path,'correlation_sd_rest')

[R1,P1] = corr(ple_task',RT_sd,'type','Pearson');
[R2,P2] = corr(ple_rest',RT_sd,'type','Pearson');

R_ple = R1;
P_ple = P1;

R_FDR1 = a_multicorrect(R1,P1,'FDR');
R_FDR2 = a_multicorrect(R2,P2,'FDR');
path = 'E:\GradCPT\visualization';
yeo_visualize(R_FDR1(2:8),path,'correlation_ple_task')
yeo_visualize(R_FDR2(2:8),path,'correlation_ple_rest')

r = corrcoef(R1,R2);

R_all = [R_sd(2:8);R_sampen(2:8);R_ple(2:8)];
P_all = [P_sd(2:8);P_sampen(2:8);P_ple(2:8)];
R_FDR_all = a_multicorrect(R_all,P_all,'FDR');


% psd analysis
psd_rest1 = psd_rest(6:42,:,:);
psd_task1 = psd_task(6:42,:,:);
for ifre = 1:size(psd_rest1,1)
    for iregion = 1:8
        x_rest = squeeze(psd_rest1(ifre,iregion,:));
        x_task = squeeze(psd_task1(ifre,iregion,:));
        [r,p] = corrcoef(x_rest,RT_sd);
        R_rest(ifre,iregion) = r(2);
        [r,p] = corrcoef(x_task,RT_sd);
        R_task(ifre,iregion) = r(2);
        P_task(ifre,iregion) = p(2);
    end
end

% calculate scale-free component

freq2 = f(6:42);


for iroi = 1:size(psd_rest1,2)
    for isubj = 1:size(psd_rest1,3)
        x1 = squeeze(psd_rest1(:,iroi,isubj));
        x2 = squeeze(psd_task1(:,iroi,isubj));
        oscillation_rest(:,iroi,isubj) = regress_osci(x1,freq2);
        oscillation_task(:,iroi,isubj) = regress_osci(x2,freq2);
    end
end

for ifre = 1:size(psd_rest1,1)
    for iregion = 1:8
        x_rest = squeeze(oscillation_rest(ifre,iregion,:));
        x_task = squeeze(oscillation_task(ifre,iregion,:));
        [r,p] = corrcoef(x_rest,RT_sd);
        R_rest(ifre,iregion) = r(2);
        [r,p] = corrcoef(x_task,RT_sd);
        R_task(ifre,iregion) = r(2);
        P_task(ifre,iregion) = p(2);
    end
end

% mediation model
[paths, stats] = mediation(sd_task(2,:)', RT_sd, ple_task(2,:)', 'plots', 'verbose', 'boot', 'bootsamples', 10000,'names', {'BOLD SD' 'RT SD' 'PLE'})
[paths, stats] = mediation(sampen_task(2,:)', RT_sd,ple_task(2,:)',  'plots', 'verbose', 'boot', 'bootsamples', 10000, 'names',{'BOLD SE' 'RT SD' 'PLE'})
[paths, stats] = mediation(sampen_rest(2,:)', RT_sd, ple_task(2,:)', 'plots', 'verbose', 'boot', 'bootsamples', 10000, 'names', {'BOLD SE' 'RT SD' 'PLE'})

[paths, stats] = mediation(sd_task(1,:)', RT_sd, ple_task(1,:)', 'plots', 'verbose', 'boot', 'bootsamples', 10000,'names', {'BOLD SD' 'RT SD' 'PLE'})
[paths, stats] = mediation(sampen_task(1,:)', RT_sd,ple_task(1,:)',  'plots', 'verbose', 'boot', 'bootsamples', 10000, 'names',{'BOLD SE' 'RT SD' 'PLE'})
[paths, stats] = mediation(sampen_rest(1,:)', RT_sd, ple_task(1,:)', 'plots', 'verbose', 'boot', 'bootsamples', 10000, 'names', {'BOLD SE' 'RT SD' 'PLE'})

% 
[paths, stats] = mediation(sd_task(2,:)', sampen_task(2,:)', ple_task(2,:)', 'plots', 'verbose', 'boot', 'bootsamples', 10000,'names', {'SD' 'SE' 'PLE'})
[paths, stats] = mediation(sd_rest(2,:)', sampen_rest(2,:)', ple_rest(1,:)', 'plots', 'verbose', 'boot', 'bootsamples', 10000,'names', {'SD' 'SE' 'PLE'})

[paths, stats] = mediation(sampen_task(2,:)',sd_task(2,:)', ple_task(2,:)', 'plots', 'verbose', 'boot', 'bootsamples', 10000,'names', {'SE' 'SD' 'PLE'})
[paths, stats] = mediation(sampen_rest(2,:)', sd_rest(2,:)',ple_rest(2,:)', 'plots', 'verbose', 'boot', 'bootsamples', 10000,'names', {'SE' 'SD' 'PLE'})

% supplement
[paths, stats] = mediation(sampen_task(2,:)', ple_task(2,:)',sd_task(2,:)', 'plots', 'verbose', 'boot', 'bootsamples', 10000,'names', {'SE' 'SD' 'PLE'})
[paths, stats] = mediation(sampen_rest(2,:)', ple_rest(2,:)',sd_rest(2,:)', 'plots', 'verbose', 'boot', 'bootsamples', 10000,'names', {'SE' 'SD' 'PLE'})

[paths, stats] = mediation(ple_task(2,:)',sampen_task(2,:)',sd_task(2,:)', 'plots', 'verbose', 'boot', 'bootsamples', 10000,'names', {'SE' 'SD' 'PLE'})
[paths, stats] = mediation(ple_rest(2,:)',sampen_rest(2,:)', sd_rest(2,:)', 'plots', 'verbose', 'boot', 'bootsamples', 10000,'names', {'SE' 'SD' 'PLE'})

[paths, stats] = mediation(ple_task(2,:)',sd_task(2,:)', sampen_task(2,:)', 'plots', 'verbose', 'boot', 'bootsamples', 10000,'names', {'SE' 'SD' 'PLE'})
[paths, stats] = mediation(ple_rest(2,:)',sd_rest(2,:)', sampen_rest(2,:)', 'plots', 'verbose', 'boot', 'bootsamples', 10000,'names', {'SE' 'SD' 'PLE'})

