%% Figure 1: Methodological validation
% gender detection
path = '/Users/yujia/Data/Data/HCP/unrestricted_ao8223730_5_24_2023_2_1_16.csv';
sub_info = readmatrix(path);



clear

path = '/Users/yujia/Data/Data/HCP/Results/mean_phase_3T';
sub_file = dir(path);
sub_file(1:2) = [];

template_path1 = '/Users/yujia/Data/Code/Common/Mask/Glasser template/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
template_path2 = '/Users/yujia/Data/Code/Common/Mask/Glasser template/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
ciftipath = '/Users/yujia/Data/Code/Common/Mask/BN_Atlas_freesurfer/fsaverage/fsaverage_LR32k/fsaverage.BN_Atlas.32k_fs_LR.dlabel.nii';
glasser_L1 = ciftiopen(template_path1);
glasser_R1 = ciftiopen(template_path2);
cifti1 = ciftiopen(ciftipath);

% compose all the subjects data to one variable
for isub = 1:length(sub_file)
    load(fullfile(path,sub_file(isub).name));
    phase_mean_stat2(:,:,isub) = mean_phase_mean_stat;
    phase_mean_stephan2(:,:,isub) =  1./mean_phase_mean_stephan;
    phase_mean_stat_shuffled2(:,:,isub) = mean_phase_mean_stat_shuffled;
end

% zscore data
for isub = 1:length(sub_file)
    for iband = 1:8

        data = squeeze(phase_mean_stat2(:,iband,isub));
        phase_mean_z(:,iband,isub) = zscore(data);

        data = squeeze(phase_mean_stephan2(:,iband,isub));
        phase_mean2_z(:,iband,isub) = zscore(data);

        data = squeeze(phase_mean_stat_shuffled2(:,iband,isub));
        phase_mean3_z(:,iband,isub) = zscore(data);
    end
end

% calculate the mean of all subjects
mean_phase_mean_z = mean(phase_mean_z,3);
mean_phase_mean2_z = mean(phase_mean2_z,3);
mean_phase_mean3_z = mean(phase_mean3_z,3);

% spatial correlation
for iband = 1:8
        cohen = squeeze(mean_phase_mean_z(:,iband));
        stephan = squeeze(mean_phase_mean2_z(:,iband));
        shuffle = squeeze(mean_phase_mean3_z(:,iband));
        [r,p] = corrcoef(cohen,stephan);
        R_mean(iband) = r(2);

        [r,p] = corrcoef(cohen,shuffle);
        R_shuffle(iband) = r(2);

end

    

for iband = 1:8
    for iroi = 1:360
        cohen = squeeze(phase_mean_z(iroi,iband,:));
        stephan = squeeze(phase_mean2_z(iroi,iband,:));
        shuffle = squeeze(phase_mean3_z(iroi,iband,:));
        icc_mean(iroi,iband) = ICC([cohen,stephan],'A-k');
        icc_shuffle(iroi,iband) = ICC([cohen,shuffle],'A-k');
    end
end

mean_icc = mean(icc_shuffle,1);

% visualize cohen and stephan result
mean_phase_mean_stat = mean(phase_mean_stat2,3);
mean_phase_mean_stephan = mean(phase_mean_stephan2,3);
mean_phase_mean_shuffled = mean(phase_mean_stat_shuffled2,3);
path = '/Users/yujia/Data/Data/HCP/Results/Visualization';


for iband = 1:2:8

    name = ['phase_mean_stat_band',num2str(iband)];
    a_visualphasevar(mean_phase_mean_stat(:,iband),cifti1,glasser_L1,glasser_R1,path,name);
    name = ['phase_mean_stephan_band',num2str(iband)];
    a_visualphasevar(mean_phase_mean_stephan(:,iband),cifti1,glasser_L1,glasser_R1,path,name);
%     name = ['phase_mean_shuffled_band',num2str(iband)];
%     a_visualphasevar(mean_phase_mean_shuffled(:,iband),cifti1,glasser_L1,glasser_R1,path,name);
end



%% Figure 2: ANOVA for 4 states

clc
clear

template_path1 = '/Users/yujia/Data/Code/Common/Mask/Glasser template/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
template_path2 = '/Users/yujia/Data/Code/Common/Mask/Glasser template/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
ciftipath = '/Users/yujia/Data/Code/Common/Mask/BN_Atlas_freesurfer/fsaverage/fsaverage_LR32k/fsaverage.BN_Atlas.32k_fs_LR.dlabel.nii';
glasser_L1 = ciftiopen(template_path1);
glasser_R1 = ciftiopen(template_path2);
cifti1 = ciftiopen(ciftipath);

path = '/Users/yujia/Data/Data/HCP/Results/mean_phase_3T';
sub_file = dir(path);
sub_file(1:2) = [];

% compose all the subjects data to one variable
for isub = 1:length(sub_file)
    load(fullfile(path,sub_file(isub).name));
    phase_mean_dy2(:,:,:,isub) = mean_phase_mean_dy;
end

phase_mean_trough2 = squeeze(phase_mean_dy2(:,:,1,:));
phase_mean_rise2 = squeeze(phase_mean_dy2(:,:,90,:));
phase_mean_peak2 = squeeze(phase_mean_dy2(:,:,180,:));
phase_mean_fall2 = squeeze(phase_mean_dy2(:,:,270,:));

trough = mean(phase_mean_trough2(:,[1,3,5,7],:),3);

rise = mean(phase_mean_rise2(:,[1,3,5,7],:),3);

peak = mean(phase_mean_peak2(:,[1,3,5,7],:),3);

fall = mean(phase_mean_fall2(:,[1,3,5,7],:),3);

all_phase = [peak,trough,rise,fall];



fre1 = [peak(:,1),trough(:,1),rise(:,1),fall(:,1)];

fre2 = [peak(:,2),trough(:,2),rise(:,2),fall(:,2)];

fre3 = [peak(:,3),trough(:,3),rise(:,3),fall(:,3)];

fre4 = [peak(:,4),trough(:,4),rise(:,4),fall(:,4)];
all_fre = [fre1,fre2,fre3,fre4];

correlation_map = corr(all_fre,all_fre,'type','Pearson');
cd('/Users/yujia/Data/Data/HCP/Figures/For seaborn')
save correlation_mat correlation_map



path = '/Users/yujia/Data/Data/HCP/Results/Visualization';
for iband = 1:2:8

    name = ['phase_mean_trough_band',num2str(iband)];
    a_visualphasevar(mean(phase_mean_trough2(:,iband,:),3),cifti1,glasser_L1,glasser_R1,path,name);

     name = ['phase_mean_rise_band',num2str(iband)];
    a_visualphasevar(mean(phase_mean_rise2(:,iband,:),3),cifti1,glasser_L1,glasser_R1,path,name);

     name = ['phase_mean_peak_band',num2str(iband)];
    a_visualphasevar(mean(phase_mean_peak2(:,iband,:),3),cifti1,glasser_L1,glasser_R1,path,name);

     name = ['phase_mean_fall_band',num2str(iband)];
    a_visualphasevar(mean(phase_mean_fall2(:,iband,:),3),cifti1,glasser_L1,glasser_R1,path,name);

end

for iroi = 1:360
  t = table(squeeze(phase_mean_trough2(iroi,1,:)),squeeze(phase_mean_rise2(iroi,1,:)),...
            squeeze(phase_mean_peak2(iroi,1,:)),squeeze(phase_mean_fall2(iroi,1,:)),... 
            squeeze(phase_mean_trough2(iroi,3,:)),squeeze(phase_mean_rise2(iroi,3,:)),...
            squeeze(phase_mean_peak2(iroi,3,:)),squeeze(phase_mean_fall2(iroi,3,:)),...  
            squeeze(phase_mean_trough2(iroi,5,:)),squeeze(phase_mean_rise2(iroi,5,:)),...
            squeeze(phase_mean_peak2(iroi,5,:)),squeeze(phase_mean_fall2(iroi,5,:)),...  
            squeeze(phase_mean_trough2(iroi,7,:)),squeeze(phase_mean_rise2(iroi,7,:)),...
            squeeze(phase_mean_peak2(iroi,7,:)),squeeze(phase_mean_fall2(iroi,7,:)),...            
            'VariableNames',...
            {'y1','y2','y3','y4','y5','y6','y7','y8','y9','y10','y11','y12',...
            'y13','y14','y15','y16'});
    meas = table({'T','R','P','F','T','R','P','F','T','R','P','F','T','R','P','F'}',...
        {'0.02 Hz','0.02 Hz','0.02 Hz','0.02 Hz',...
        '0.04 Hz','0.04 Hz','0.04 Hz','0.04 Hz',...
        '0.06 Hz','0.06 Hz','0.06 Hz','0.06 Hz',...
        '0.08 Hz','0.08 Hz','0.08 Hz','0.08 Hz'}',...
        'VariableNames',{'phase','frequency'});
    rm = fitrm(t,'y1-y16~1','WithinDesign',meas);
    tbl_ep = epsilon(rm);
    [tbl,A,C,D] = ranova(rm,'WithinModel',...
        'phase+frequency+phase*frequency');
    F_phase(iroi) = tbl{3,4};
    F_frequency(iroi) = tbl{5,4};
    F_inter(iroi) = tbl{7,4};
    P_phase(iroi) = tbl{3,5};  
    P_frequency(iroi) = tbl{5,5}; 
    P_inter(iroi) = tbl{7,5}; 


 end

F_phase_FDR = a_multicorrect(F_phase,P_phase,'FDR');
F_frequency_FDR = a_multicorrect(F_frequency,P_frequency,'FDR');
F_inter_FDR = a_multicorrect(F_inter,P_inter,'FDR');


path = '/Users/yujia/Data/Data/HCP/Results/Visualization';
a_visualphasevar(F_phase_FDR,cifti1,glasser_L1,glasser_R1,path,'ANOVA_phase_mean');
a_visualphasevar(F_frequency_FDR,cifti1,glasser_L1,glasser_R1,path,'ANOVA_frequency_mean');
a_visualphasevar(F_inter_FDR,cifti1,glasser_L1,glasser_R1,path,'ANOVA_inter_mean');

path = '/Users/yujia/Data/Data/HCP/Results/Visualization/t-test_mean';
for iband = 1:2:8

    [h,p,ci,stats] = ttest(squeeze(phase_mean_trough2(:,iband,:))',...
        squeeze(phase_mean_rise2(:,iband,:))');
    cohend = a_cohend(stats.tstat,73);
    name = ['cohend-T-R-band',num2str(iband)];
    a_visualphasevar(cohend,cifti1,glasser_L1,glasser_R1,path,name);
    T_FDR = a_multicorrect(stats.tstat,p,'FDR');
    name = ['T-R-band',num2str(iband)];
    a_visualphasevar(T_FDR,cifti1,glasser_L1,glasser_R1,path,name);

    [h,p,ci,stats] = ttest(squeeze(phase_mean_trough2(:,iband,:))',...
        squeeze(phase_mean_peak2(:,iband,:))');
    cohend = a_cohend(stats.tstat,73);
    name = ['cohend-T-P-band',num2str(iband)];
    a_visualphasevar(cohend,cifti1,glasser_L1,glasser_R1,path,name);
    T_FDR = a_multicorrect(stats.tstat,p,'FDR');
    name = ['T-P-band',num2str(iband)];
    a_visualphasevar(T_FDR,cifti1,glasser_L1,glasser_R1,path,name);

    [h,p,ci,stats] = ttest(squeeze(phase_mean_trough2(:,iband,:))',...
        squeeze(phase_mean_fall2(:,iband,:))');
    cohend = a_cohend(stats.tstat,73);
    name = ['cohend-T-F-band',num2str(iband)];
    a_visualphasevar(cohend,cifti1,glasser_L1,glasser_R1,path,name);
    T_FDR = a_multicorrect(stats.tstat,p,'FDR');
    name = ['T-F-band',num2str(iband)];
    a_visualphasevar(T_FDR,cifti1,glasser_L1,glasser_R1,path,name);

    [h,p,ci,stats] = ttest(squeeze(phase_mean_rise2(:,iband,:))',...
        squeeze(phase_mean_peak2(:,iband,:))');
    cohend = a_cohend(stats.tstat,73);
    name = ['cohend-R-P-band',num2str(iband)];
    a_visualphasevar(cohend,cifti1,glasser_L1,glasser_R1,path,name);
    T_FDR = a_multicorrect(stats.tstat,p,'FDR');
    name = ['R-P-band',num2str(iband)];
    a_visualphasevar(T_FDR,cifti1,glasser_L1,glasser_R1,path,name);

    
    [h,p,ci,stats] = ttest(squeeze(phase_mean_rise2(:,iband,:))',...
        squeeze(phase_mean_fall2(:,iband,:))');
    cohend = a_cohend(stats.tstat,73);
    name = ['cohend-R-F-band',num2str(iband)];
    a_visualphasevar(cohend,cifti1,glasser_L1,glasser_R1,path,name);
    T_FDR = a_multicorrect(stats.tstat,p,'FDR');
    name = ['R-F-band',num2str(iband)];
    a_visualphasevar(T_FDR,cifti1,glasser_L1,glasser_R1,path,name);

    [h,p,ci,stats] = ttest(squeeze(phase_mean_peak2(:,iband,:))',...
        squeeze(phase_mean_fall2(:,iband,:))');
    cohend = a_cohend(stats.tstat,73);
    name = ['cohend-P-F-band',num2str(iband)];
    a_visualphasevar(cohend,cifti1,glasser_L1,glasser_R1,path,name);
    T_FDR = a_multicorrect(stats.tstat,p,'FDR');
    name = ['P-F-band',num2str(iband)];
    a_visualphasevar(T_FDR,cifti1,glasser_L1,glasser_R1,path,name);
end

% for python
band=[1,3,5,7];


x1 = squeeze(mean(phase_mean_fall2(:,band,:),1));
x_fall = squeeze(mean(x1,1));
x_fall_allband = x1;
x1 = squeeze(mean(phase_mean_peak2(:,band,:),1));
x_peak = squeeze(mean(x1,1));
x_peak_allband = x1;
x1 = squeeze(mean(phase_mean_rise2(:,band,:),1));
x_rise = squeeze(mean(x1,1));
x_rise_allband = x1;
x1 = squeeze(mean(phase_mean_trough2(:,band,:),1));
x_trough = squeeze(mean(x1,1));
x_trough_allband = x1;



t = table(squeeze(x_trough_allband(1,:))',squeeze(x_rise_allband(1,:))',...
        squeeze(x_peak_allband(1,:))',squeeze(x_fall_allband(1,:))',... 
        squeeze(x_trough_allband(2,:))',squeeze(x_rise_allband(2,:))',...
        squeeze(x_peak_allband(2,:))',squeeze(x_fall_allband(2,:))',...
        squeeze(x_trough_allband(3,:))',squeeze(x_rise_allband(3,:))',...
        squeeze(x_peak_allband(3,:))',squeeze(x_fall_allband(3,:))',...
        squeeze(x_trough_allband(4,:))',squeeze(x_rise_allband(4,:))',...
        squeeze(x_peak_allband(4,:))',squeeze(x_fall_allband(4,:))',...
        'VariableNames',...
        {'y1','y2','y3','y4','y5','y6','y7','y8','y9','y10','y11','y12',...
        'y13','y14','y15','y16'});
meas = table({'T','R','P','F','T','R','P','F','T','R','P','F','T','R','P','F'}',...
    {'0.02 Hz','0.02 Hz','0.02 Hz','0.02 Hz',...
    '0.04 Hz','0.04 Hz','0.04 Hz','0.04 Hz',...
    '0.06 Hz','0.06 Hz','0.06 Hz','0.06 Hz',...
    '0.08 Hz','0.08 Hz','0.08 Hz','0.08 Hz'}',...
    'VariableNames',{'phase','frequency'});
rm = fitrm(t,'y1-y16~1','WithinDesign',meas);
tbl_ep = epsilon(rm);
[tbl,A,C,D] = ranova(rm,'WithinModel',...
    'phase+frequency+phase*frequency'); 

x_band1 = [x_peak_allband(1,:),x_trough_allband(1,:),x_rise_allband(1,:),x_fall_allband(1,:)]';
x_band2 = [x_peak_allband(2,:),x_trough_allband(2,:),x_rise_allband(2,:),x_fall_allband(2,:)]';
x_band3 = [x_peak_allband(3,:),x_trough_allband(3,:),x_rise_allband(3,:),x_fall_allband(3,:)]';
x_band4 = [x_peak_allband(4,:),x_trough_allband(4,:),x_rise_allband(4,:),x_fall_allband(4,:)]';


for iband = 1:4
    y=[x_peak_allband(iband,:)',x_trough_allband(iband,:)',x_rise_allband(iband,:)',x_fall_allband(iband,:)'];
    [p,tbl,stats] = anova1(y,[],'off');
    F_phase(iband) = tbl{2,5};
    eta(iband) = a_etasq(F_phase(iband),3,720);
end

x_mean_phase = [x_peak,x_trough,x_rise,x_fall]';


x2 = cat(4,phase_mean_fall2(:,band,:),phase_mean_peak2(:,band,:),...
    phase_mean_rise2(:,band,:),phase_mean_trough2(:,band,:));

x3 = squeeze(mean(x2,1));
x_mean_frequency = squeeze(mean(x3,3))';

x_mean_frequency = reshape(x_mean_frequency,[],1);



cd('/Users/yujia/Data/Data/HCP/Figures/For seaborn')
save figure2_bar x_mean_phase x_mean_frequency  x_band1 x_band2 x_band3 x_band4




%% Figure 3: Phase dynamics and modulation index
clc
clear

template_path1 = 'D:\Code\Common\Mask\Glasser template\Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
template_path2 = 'D:\Code\Common\Mask\Glasser template\Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
ciftipath = 'D:\Code\Common\Mask\BN_Atlas_freesurfer\fsaverage\fsaverage_LR32k\fsaverage.BN_Atlas.32k_fs_LR.dlabel.nii';
glasser_L1 = ciftiopen(template_path1);
glasser_R1 = ciftiopen(template_path2);
cifti1 = ciftiopen(ciftipath);

path = 'E:\Data\HCP\Results\mean_phase_3T';
sub_file = dir(path);   
sub_file(1:2) = [];

for isub = 1:length(sub_file)
    load(fullfile(path,sub_file(isub).name));
    phase_mean_dy2(:,:,:,isub) = mean_phase_mean_dy;
end

% line chart
mean_phase_mean_dy = squeeze(mean(phase_mean_dy2,1));

mean_phase_dy = mean(mean_phase_mean_dy,3);

for iband = 1:2:8
    KL_all(iband) = a_mi(mean_phase_dy(iband,:));
end



phase_num = -180:179;
phase_num = phase_num';


% for python 
x = permute(mean_phase_mean_dy,[2 3 1]);
x = reshape(x ,[],8);
x = x(:,[1,3,5,7]);


mean_x = squeeze(mean(mean_phase_mean_dy,3))';
KL_mean = a_mi(mean_x(:,1));


y = repmat(phase_num,181,1);

cd('/Users/yujia/Data/Data/HCP/Figures/For seaborn')
save phase_mean_dy x y


% KL all roi bar

mean_phase_mean_dy = squeeze(mean(phase_mean_dy2,1));

for isub = 1:size(mean_phase_mean_dy,3)
    for iband = 1:8
        data = mean_phase_mean_dy(iband,:,isub)';
        score_KL_mean(iband,isub) = a_mi(data);
    end
end

t = table(squeeze(score_KL_mean(1,:))',squeeze(score_KL_mean(3,:))',...
        squeeze(score_KL_mean(5,:))',squeeze(score_KL_mean(7,:))',...            
            'VariableNames',{'y1','y2','y3','y4'});
meas = table({'0.02','0.04','0.06','0.08'}',... 
    'VariableNames',{'frequency'});
rm = fitrm(t,'y1-y4~1','WithinDesign',meas);
tbl_ep = epsilon(rm);
[tbl,A,C,D] = ranova(rm,'WithinModel',...
    'frequency');

% for python
x2 = score_KL_mean([1,3,5,7],:)';

x2 = reshape(x2,[],1);

cd('/Users/yujia/Data/Data/HCP/Figures/For seaborn')
save figure3_bar x2 

% topography and anova
path = '/Users/yujia/Data/Data/HCP/Results/mean_phase_3T';
sub_file = dir(path);
sub_file(1:2) = [];

for isub = 1:length(sub_file)
    isub
    load(fullfile(path,sub_file(isub).name));
    for iroi = 1:360
        for iband = 1:8
            data = squeeze(mean_phase_mean_dy(iroi,iband,:));
            score_KL_mean_topo(iroi,iband,isub) = a_mi(data);
        end
    end
end

mean_mi_mean = mean(score_KL_mean_topo,3);

path = '/Users/yujia/Data/Data/HCP/Results/Visualization';
for iband = 1:2:8
   
    name = ['MI_mean_band',num2str(iband)];
    a_visualphasevar(mean_mi_mean(:,iband),cifti1,glasser_L1,glasser_R1,path,name);
end

% anova
for iroi = 1:360

    t = table(squeeze(score_KL_mean_topo(iroi,1,:)),squeeze(score_KL_mean_topo(iroi,3,:)),...
        squeeze(score_KL_mean_topo(iroi,5,:)),squeeze(score_KL_mean_topo(iroi,7,:)),...            
            'VariableNames',{'y1','y2','y3','y4'});
    meas = table({'0.02','0.04','0.06','0.08'}',... 
        'VariableNames',{'frequency'});
    rm = fitrm(t,'y1-y4~1','WithinDesign',meas);
    tbl_ep = epsilon(rm);
    [tbl,A,C,D] = ranova(rm,'WithinModel',...
        'frequency');
    F(iroi) = tbl{3,4};

    P(iroi) = tbl{3,5}; 
end

F_FDR = a_multicorrect(F,P,'FDR');
path = '/Users/yujia/Data/Data/HCP/Results/Visualization';
a_visualphasevar(F_FDR,cifti1,glasser_L1,glasser_R1,path,'ANOVA_MI_mean');

%% Figure 4: ACW and frequency
clear
path = 'E:\Data\HCP\Results2\mean_acw_3T';
sub_file = dir(path);
sub_file(1:2) = [];
  
for isub = 1:length(sub_file)
    load(fullfile(path,sub_file(isub).name));
    acw_0_all(:,isub) = mean_acw;
end

mean_acw_0 = mean(acw_0_all,2);

template_path1 = '/Users/yujia/Data/Code/Common/Mask/Glasser template/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
template_path2 = '/Users/yujia/Data/Code/Common/Mask/Glasser template/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
ciftipath = '/Users/yujia/Data/Code/Common/Mask/BN_Atlas_freesurfer/fsaverage/fsaverage_LR32k/fsaverage.BN_Atlas.32k_fs_LR.dlabel.nii';
glasser_L1 = ciftiopen(template_path1);
glasser_R1 = ciftiopen(template_path2);
cifti1 = ciftiopen(ciftipath);
path = '/Users/yujia/Data/Data/HCP/Results/Visualization';


name = ['acw_0'];
a_visualphasevar(mean_acw_0,cifti1,glasser_L1,glasser_R1,path,name);

% scatter chart for mean ACW and freuqnecy

path = 'E:\Data\HCP\Results2/mean_phase_3T';
sub_file = dir(path);
sub_file(1:2) = [];

% static phase var
for isub = 1:length(sub_file)
    load(fullfile(path,sub_file(isub).name));
    phase_mean_stat2(:,:,isub) = mean_phase_mean_stat;
end

mean_phase_mean_stat2 = mean(phase_mean_stat2,3);

for iband = 1:4
    [r,p] = corrcoef(mean_acw_0,mean_phase_mean_stat2(:,iband));
    R_mean(iband) = r(2);
    P_mean(iband) = p(2);
end

% for python
x2 = mean_phase_mean_stat2;
y1 = mean_acw_0;

cd('/Users/yujia/Data/Data/HCP/Figures/For seaborn')
save figure4_scatter x2 y1 

x1 = R_mean([1,3,5,7]);
save figure4_bar x1

for i1 = 1:4
    for i2 = 1:4
        [Zscore(i1,i2), Pvalue(i1,i2)] = gretna_ztest_two_corrcoef(x1(i1), x1(i2), 360, 360,'both');
    end
end

%% compare shulffed data
clear
path = 'E:\Data\HCP\Results\mean_acw_3T';
sub_file = dir(path);
sub_file(1:2) = [];
  
for isub = 1:length(sub_file)
    load(fullfile(path,sub_file(isub).name));
    acw_0_all(:,isub) = mean_acw;
end
mean_acw_0 = mean(acw_0_all,2);


path = 'E:\Data\HCP\Results\mean_phase_3T';
sub_file = dir(path);
sub_file(1:2) = [];

% static phase var
for isub = 1:length(sub_file)
    load(fullfile(path,sub_file(isub).name));
    phase_mean_stat2(:,:,isub) = mean_phase_mean_stat_shuffled;
end

mean_phase_mean_stat2 = mean(phase_mean_stat2,3);

for iband = 1:8
    [r,p] = corrcoef(mean_acw_0,mean_phase_mean_stat2(:,iband));
    R_mean(iband) = r(2);
    P_mean(iband) = p(2);
end

%% Figure 5 task

% INT map for task and correlation with frequency sliding
clear
clc

path_rest = '/Users/yujia/Data/Data/HCP/Results/mean_acw_7T_rest';
path_task = '/Users/yujia/Data/Data/HCP/Results/mean_acw_7T_task';
sub_file = dir(path_rest);
sub_file(1:2) = [];

for isub = 1:length(sub_file)
    load(fullfile(path_rest,sub_file(isub).name));
    acw_0_rest(:,isub) = mean_acw;

    load(fullfile(path_task,sub_file(isub).name));
    acw_0_task(:,isub) = mean_acw;
end


acw_0_rest_topo = squeeze(mean(acw_0_rest,2));
acw_0_task_topo = squeeze(mean(acw_0_task,2));


template_path1 = '/Users/yujia/Data/Code/Common/Mask/Glasser template/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
template_path2 = '/Users/yujia/Data/Code/Common/Mask/Glasser template/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
ciftipath = '/Users/yujia/Data/Code/Common/Mask/BN_Atlas_freesurfer/fsaverage/fsaverage_LR32k/fsaverage.BN_Atlas.32k_fs_LR.dlabel.nii';
glasser_L1 = ciftiopen(template_path1);
glasser_R1 = ciftiopen(template_path2);
cifti1 = ciftiopen(ciftipath);
path = '/Users/yujia/Data/Data/HCP/Results/Visualization';

name = ['acw_0_rest_topo'];
a_visualphasevar(acw_0_rest_topo(:,1),cifti1,glasser_L1,glasser_R1,path,name);

name = ['acw_0_task_topo'];
a_visualphasevar(acw_0_task_topo(:,1),cifti1,glasser_L1,glasser_R1,path,name);

% static phase var
path_rest = '/Users/yujia/Data/Data/HCP/Results/mean_phase_7T_rest';
path_task = '/Users/yujia/Data/Data/HCP/Results/mean_phase_7T_task';
sub_file = dir(path_rest);
sub_file(1:2) = [];

% static phase var
for isub = 1:length(sub_file)
    load(fullfile(path_rest,sub_file(isub).name));
    phase_mean_rest_stat2(:,:,isub) = mean_phase_mean_stat;

    load(fullfile(path_task,sub_file(isub).name));
    phase_mean_task_stat2(:,:,isub) = mean_phase_mean_stat;
end


mean_phase_mean_rest_stat2 = mean(phase_mean_rest_stat2,3);
mean_phase_mean_task_stat2 = mean(phase_mean_task_stat2,3);

for iband = 1:8
    [r,p] = corrcoef(mean_phase_mean_rest_stat2(:,iband),acw_0_rest_topo);
    R_mean_rest(iband) = r(2);
    [r,p] = corrcoef(mean_phase_mean_task_stat2(:,iband),acw_0_task_topo);
    R_mean_task(iband) = r(2);
end

template_path1 = '/Users/yujia/Data/Code/Common/Mask/Glasser template/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
template_path2 = '/Users/yujia/Data/Code/Common/Mask/Glasser template/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
ciftipath = '/Users/yujia/Data/Code/Common/Mask/BN_Atlas_freesurfer/fsaverage/fsaverage_LR32k/fsaverage.BN_Atlas.32k_fs_LR.dlabel.nii';
glasser_L1 = ciftiopen(template_path1);
glasser_R1 = ciftiopen(template_path2);
cifti1 = ciftiopen(ciftipath);
path = '/Users/yujia/Data/Data/HCP/Results/Visualization';
for iband = 1:8    
    name = ['mean_phase_mean_rest_band',num2str(iband)];
    a_visualphasevar(mean_phase_mean_rest_stat2(:,1),cifti1,glasser_L1,glasser_R1,path,name);
    name = ['mean_phase_mean_task_band',num2str(iband)];
    a_visualphasevar(mean_phase_mean_rest_stat2(:,1),cifti1,glasser_L1,glasser_R1,path,name);
end

%% Figure 6 rest-task
clear
clc

path_rest = '/Users/yujia/Data/Data/HCP/Results/mean_phase_7T_rest';
path_task = '/Users/yujia/Data/Data/HCP/Results/mean_phase_7T_task';
sub_file = dir(path_rest);
sub_file(1:2) = [];

% compose all the subjects data to one variable
for isub = 1:length(sub_file)
    load(fullfile(path_rest,sub_file(isub).name));
    phase_mean_rest_dy2(:,:,:,isub) = mean_phase_mean_dy;

    load(fullfile(path_task,sub_file(isub).name));
    phase_mean_task_dy2(:,:,:,isub) = mean_phase_mean_dy;
end

phase_trough_rest = squeeze(phase_mean_rest_dy2(:,:,1,:));
phase_peak_rest = squeeze(phase_mean_rest_dy2(:,:,180,:));
phase_rise_rest = squeeze(phase_mean_rest_dy2(:,:,90,:));
phase_fall_rest = squeeze(phase_mean_rest_dy2(:,:,270,:));

phase_trough_task = squeeze(phase_mean_task_dy2(:,:,1,:));
phase_peak_task = squeeze(phase_mean_task_dy2(:,:,180,:));
phase_rise_task = squeeze(phase_mean_task_dy2(:,:,90,:));
phase_fall_task = squeeze(phase_mean_task_dy2(:,:,270,:));

% for python
cd('/Users/yujia/Data/Data/HCP/Figures/For seaborn')
mean_phase_trough_rest = squeeze(mean(phase_trough_rest(:,1,:),1));
mean_phase_trough_task = squeeze(mean(phase_trough_task(:,1,:),1));
mean_phase_peak_rest = squeeze(mean(phase_peak_rest(:,1,:),1));
mean_phase_peak_task = squeeze(mean(phase_peak_task(:,1,:),1));
mean_phase_rise_rest = squeeze(mean(phase_rise_rest(:,1,:),1));
mean_phase_rise_task = squeeze(mean(phase_rise_task(:,1,:),1));
mean_phase_fall_rest = squeeze(mean(phase_fall_rest(:,1,:),1));
mean_phase_fall_task = squeeze(mean(phase_fall_task(:,1,:),1));


x = [mean_phase_peak_rest;
    mean_phase_peak_task;mean_phase_trough_rest;mean_phase_trough_task;mean_phase_rise_rest;mean_phase_rise_task;
    mean_phase_fall_rest;mean_phase_fall_task];

save figure6_fs x

% t-test for fs
[h,p,ci,stats] = ttest(mean_phase_peak_rest,mean_phase_peak_task);
d1 = a_cohend(stats.tstat,179);
[h,p,ci,stats] = ttest(mean_phase_trough_rest,mean_phase_trough_task);
d2 = a_cohend(stats.tstat,179);
[h,p,ci,stats] = ttest(mean_phase_rise_rest,mean_phase_rise_task);
d3 = a_cohend(stats.tstat,179);
[h,p,ci,stats] = ttest(mean_phase_fall_rest,mean_phase_fall_task);
d4 = a_cohend(stats.tstat,179);


for iband = 1:8
    for iroi = 1:360

        trough_rest = squeeze(phase_trough_rest(iroi,iband,:));
        trough_task = squeeze(phase_trough_task(iroi,iband,:));
        peak_rest = squeeze(phase_peak_rest(iroi,iband,:));
        peak_task = squeeze(phase_peak_task(iroi,iband,:));
        rise_rest = squeeze(phase_rise_rest(iroi,iband,:));
        rise_task = squeeze(phase_rise_task(iroi,iband,:));
        fall_rest = squeeze(phase_fall_rest(iroi,iband,:));
        fall_task = squeeze(phase_fall_task(iroi,iband,:));

        [h,p,ci,stats] = ttest(trough_rest,trough_task);
        P_trough(iroi,iband) = p;
        T_trough(iroi,iband) = stats.tstat;
        [h,p,ci,stats] = ttest(peak_rest,peak_task);
        P_peak(iroi,iband) = p;
        T_peak(iroi,iband) = stats.tstat;
        [h,p,ci,stats] = ttest(rise_rest,rise_task);
        P_rise(iroi,iband) = p;
        T_rise(iroi,iband) = stats.tstat;
        [h,p,ci,stats] = ttest(fall_rest,fall_task);
        P_fall(iroi,iband) = p;
        T_fall(iroi,iband) = stats.tstat;


    end    
end

for iband = 1:8
    T_trough_correct(:,iband) = a_multicorrect(T_trough(:,iband),P_trough(:,iband),'FDR');
    T_peak_correct(:,iband) = a_multicorrect(T_peak(:,iband),P_peak(:,iband),'FDR');
    T_rise_correct(:,iband) = a_multicorrect(T_rise(:,iband),P_rise(:,iband),'FDR');
    T_fall_correct(:,iband) = a_multicorrect(T_fall(:,iband),P_fall(:,iband),'FDR');
end

template_path1 = '/Users/yujia/Data/Code/Common/Mask/Glasser template/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
template_path2 = '/Users/yujia/Data/Code/Common/Mask/Glasser template/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
ciftipath = '/Users/yujia/Data/Code/Common/Mask/BN_Atlas_freesurfer/fsaverage/fsaverage_LR32k/fsaverage.BN_Atlas.32k_fs_LR.dlabel.nii';
glasser_L1 = ciftiopen(template_path1);
glasser_R1 = ciftiopen(template_path2);
cifti1 = ciftiopen(ciftipath);
path = '/Users/yujia/Data/Data/HCP/Results/Visualization';

name = ['T_trough_rest_task_band',num2str(1)];
a_visualphasevar(T_trough_correct(:,1),cifti1,glasser_L1,glasser_R1,path,name);
name = ['T_peak_rest_task_band',num2str(1)];
a_visualphasevar(T_peak_correct(:,1),cifti1,glasser_L1,glasser_R1,path,name);
name = ['T_rise_rest_task_band',num2str(1)];
a_visualphasevar(T_rise_correct(:,1),cifti1,glasser_L1,glasser_R1,path,name);
name = ['T_fall_rest_task_band',num2str(1)];
a_visualphasevar(T_fall_correct(:,1),cifti1,glasser_L1,glasser_R1,path,name);



% KL divergence
path_rest = '/Users/yujia/Data/Data/HCP/Results/mean_phase_7T_rest';
path_task = '/Users/yujia/Data/Data/HCP/Results/mean_phase_7T_task';
sub_file = dir(path_rest);
sub_file(1:2) = [];


for isub = 1:length(sub_file)
    load(fullfile(path_rest,sub_file(isub).name));
    phase_mean_dy2_rest(:,:,:,isub) = mean_phase_mean_dy;

    load(fullfile(path_task,sub_file(isub).name));
    phase_mean_dy2_task(:,:,:,isub) = mean_phase_mean_dy;
end

mean_phase_mean_dy_rest = squeeze(mean(phase_mean_dy2_rest,1));
mean_phase_mean_dy_task = squeeze(mean(phase_mean_dy2_task,1));

% 所有roi平均
for isub = 1:size(mean_phase_mean_dy_rest,3)
    for iband = 1:8
        data = mean_phase_mean_dy_rest(iband,:,isub)';
        score_KL_mean_rest_all(iband,isub) = a_mi(data);

        data = mean_phase_mean_dy_task(iband,:,isub)';
        score_KL_mean_task_all(iband,isub) = a_mi(data);

    end
end

% for python
A = score_KL_mean_rest_all;
B = score_KL_mean_task_all;
A([2,4,6,8],:) = [];
B([2,4,6,8],:) = [];
x1 = reshape(A',[],1);
x2 = reshape(B',[],1);
x = [x1;x2];
cd('/Users/yujia/Data/Data/HCP/Figures/For seaborn')
save figure6_kl x

[h,p,ci,stats] = ttest(score_KL_mean_rest_all(1,:)',score_KL_mean_task_all(1,:)');


for iband = 1:8

        mean_rest = squeeze(score_KL_mean_rest_all(iband,:));
        mean_task = squeeze(score_KL_mean_task_all(iband,:));

        [h,p,ci,stats] = ttest2(mean_rest,mean_task);
        P_mean_all(iband) = p;
        T_mean_all(iband) = stats.tstat;
end


% 拓扑图
for isub = 1:length(sub_file)
    isub
    load(fullfile(path_rest,sub_file(isub).name));
    for iroi = 1:360
        for iband = 1:8
            data = squeeze(mean_phase_mean_dy(iroi,iband,:));
            score_KL_mean_rest(iroi,iband,isub) = a_mi(data);
        end
    end
end


for isub = 1:length(sub_file)
    isub
    load(fullfile(path_task,sub_file(isub).name));
    for iroi = 1:360
        for iband = 1:8
            data = squeeze(mean_phase_mean_dy(iroi,iband,:));
            score_KL_mean_task(iroi,iband,isub) = a_mi(data);
        end
    end
end

% t-test
for iband = 1:8
    for iroi = 1:360

        mean_rest = squeeze(score_KL_mean_rest(iroi,iband,:));
        mean_task = squeeze(score_KL_mean_task(iroi,iband,:));


        [h,p,ci,stats] = ttest(mean_rest,mean_task);
        P_mean(iroi,iband) = p;
        T_mean(iroi,iband) = stats.tstat;

    end    
end

for iband = 1:8
    T_mean_correct(:,iband) = a_multicorrect(T_mean(:,iband),P_mean(:,iband),'FDR');
  
end

template_path1 = '/Users/yujia/Data/Code/Common/Mask/Glasser template/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
template_path2 = '/Users/yujia/Data/Code/Common/Mask/Glasser template/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
ciftipath = '/Users/yujia/Data/Code/Common/Mask/BN_Atlas_freesurfer/fsaverage/fsaverage_LR32k/fsaverage.BN_Atlas.32k_fs_LR.dlabel.nii';
glasser_L1 = ciftiopen(template_path1);
glasser_R1 = ciftiopen(template_path2);
cifti1 = ciftiopen(ciftipath);
path = '/Users/yujia/Data/Data/HCP/Results/Visualization';

name = ['T_phase_mean_rest_task_band1'];
a_visualphasevar(T_mean_correct(:,1),cifti1,glasser_L1,glasser_R1,path,name);

name = ['T_phase_mean_rest_task_band1_unthre'];
a_visualphasevar(T_mean(:,1),cifti1,glasser_L1,glasser_R1,path,name);



% INT rest-task difference

path_rest = '/Users/yujia/Data/Data/HCP/Results/mean_acw_7T_rest';
path_task = '/Users/yujia/Data/Data/HCP/Results/mean_acw_7T_task';
sub_file = dir(path_rest);
sub_file(1:2) = [];

for isub = 1:length(sub_file)
    load(fullfile(path_rest,sub_file(isub).name));
    acw_0_rest(:,isub) = mean_acw;

    load(fullfile(path_task,sub_file(isub).name));
    acw_0_task(:,isub) = mean_acw;
end

mean_acw_0_rest = squeeze(mean(acw_0_rest,1));
mean_acw_0_task = squeeze(mean(acw_0_task,1));

% for python

x_mean_rest = score_KL_mean_rest(1,:)';
x_mean_task = score_KL_mean_task(1,:)';


x1 = [x_mean_rest;x_mean_task];

x3 = [mean_acw_0_rest';mean_acw_0_task'];

cd('/Users/yujia/Data/Data/HCP/Figures/For seaborn')
save figure5_bar x1  x3





[h,p,ci,stats] = ttest(mean_acw_0_rest,mean_acw_0_task);

for iroi = 1:360
    rest_0 = acw_0_rest(iroi,:);
    task_0 = acw_0_task(iroi,:);
    [h,p,ci,stats] = ttest(rest_0,task_0);
    P_0(iroi) = p;
    T_0(iroi) = stats.tstat;

end

T_0_FDR = a_multicorrect(T_0,P_0,'FDR');

template_path1 = '/Users/yujia/Data/Code/Common/Mask/Glasser template/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
template_path2 = '/Users/yujia/Data/Code/Common/Mask/Glasser template/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
ciftipath = '/Users/yujia/Data/Code/Common/Mask/BN_Atlas_freesurfer/fsaverage/fsaverage_LR32k/fsaverage.BN_Atlas.32k_fs_LR.dlabel.nii';
glasser_L1 = ciftiopen(template_path1);
glasser_R1 = ciftiopen(template_path2);
cifti1 = ciftiopen(ciftipath);
path = '/Users/yujia/Data/Data/HCP/Results/Visualization';

name = ['T_acw_rest_task_FDR'];
a_visualphasevar(T_0_FDR,cifti1,glasser_L1,glasser_R1,path,name);




path_rest = '/Users/yujia/Data/Data/HCP/Results/mean_phase_7T_rest';
path_task = '/Users/yujia/Data/Data/HCP/Results/mean_phase_7T_task';
sub_file = dir(path_rest);
sub_file(1:2) = [];

% compose all the subjects data to one variable
for isub = 1:length(sub_file)
    load(fullfile(path_rest,sub_file(isub).name));
    phase_mean_rest_stat2(:,:,isub) = mean_phase_mean_stat;

    load(fullfile(path_task,sub_file(isub).name));
    phase_mean_task_stat2(:,:,isub) = mean_phase_mean_stat;
end


for iband = 1:8
    for iroi = 1:360

        rest = squeeze(phase_mean_rest_stat2(iroi,iband,:));
        task = squeeze(phase_mean_task_stat2(iroi,iband,:));


        [h,p,ci,stats] = ttest(rest,task);
        P_phase(iroi,iband) = p;
        T_phase(iroi,iband) = stats.tstat;

    end    
end


for i = 1:8
    [r,p] = corrcoef(T_0',T_phase(:,i));
    R_rest_task_diff(i) = r(2);
    P_rest_task_diff(i) = p(2);
end

% for python
cd('/Users/yujia/Data/Data/HCP/Figures/For seaborn') 

x1 = R_rest_task_diff([1,3,5,7]);
save figure7_bar x1


path = '/Users/yujia/Data/Data/HCP/Results/Visualization';

for i = 1:8
    name = ['T_fs_stat_rest_task_band',num2str(i)];
    a_visualphasevar(T_phase(:,i),cifti1,glasser_L1,glasser_R1,path,name);
end



%% Supplement for Georg: ACW, FS in different state
clear

path = '/Users/yujia/Data/Data/HCP/Results/mean_acw_7T_task';
sub_file = dir(path);
sub_file(1:2) = [];
  
for isub = 1:length(sub_file)
    load(fullfile(path,sub_file(isub).name));
    acw_0_all(:,isub) = mean_acw;
end

mean_acw_0_task = mean(acw_0_all,2);


% scatter chart for mean ACW and freuqnecy

path = '/Users/yujia/Data/Data/HCP/Results/mean_phase_7T_task';
sub_file = dir(path);
sub_file(1:2) = [];

% static phase var
for isub = 1:length(sub_file)
    load(fullfile(path,sub_file(isub).name));
    phase_mean_stat2(:,:,isub) = mean_phase_mean_stat;
end

mean_phase_mean_stat2_task = mean(phase_mean_stat2,3);


path = '/Users/yujia/Data/Data/HCP/Results/mean_acw_7T_rest';
sub_file = dir(path);
sub_file(1:2) = [];
  
for isub = 1:length(sub_file)
    load(fullfile(path,sub_file(isub).name));
    acw_0_all(:,isub) = mean_acw;
end

mean_acw_0_rest = mean(acw_0_all,2);

% scatter chart for mean ACW and freuqnecy

path = '/Users/yujia/Data/Data/HCP/Results/mean_phase_7T_rest';
sub_file = dir(path);
sub_file(1:2) = [];

% static phase var
for isub = 1:length(sub_file)
    load(fullfile(path,sub_file(isub).name));
    phase_mean_stat2(:,:,isub) = mean_phase_mean_stat;
end

mean_phase_mean_stat2_rest = mean(phase_mean_stat2,3);

for iband = 1:8
    [r,p] = corrcoef(mean_acw_0_rest,mean_phase_mean_stat2_task(:,iband));
    R_mean(iband) = r(2);
    P(iband) = p(2);
end
    
% for python
x2 = mean_phase_mean_stat2;
y1 = mean_acw_0;

cd('/Users/yujia/Data/Data/HCP/Figures/For seaborn 7T task')
save figure4_scatter x2 y1 

x1 = R_mean([1,3,5,7]);
save figure4_bar x1

for i1 = 1:4
    for i2 = 1:4
        [Zscore(i1,i2), Pvalue(i1,i2)] = gretna_ztest_two_corrcoef(x1(i1), x1(i2), 360, 360,'both');
    end
end



% scatter chart for mean ACW and freuqnecy

path = 'E:\Phase dynamics and INT\Results\mean_fp_3T';
sub_file = dir(path);
sub_file(1:2) = [];

% static phase var
for isub = 1:length(sub_file)
    load(fullfile(path,sub_file(isub).name));
    phase_mean_stat2(:,:,isub) = mean_fp;
end

mean_phase_mean_stat2 = mean(phase_mean_stat2,3);

for iband = 1:2:8
    [r,p] = corrcoef(mean_acw_0,mean_phase_mean_stat2(:,iband));
    R_mean(iband) = r(2);
    P_mean(iband) = p(2);
end

R_mean([2,4,6]) = [];

template_path1 = 'D:\Code\Common\Mask\Glasser template\Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
template_path2 = 'D:\Code\Common\Mask\Glasser template\Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
ciftipath = 'D:\Code\Common\Mask\BN_Atlas_freesurfer\fsaverage\fsaverage_LR32k\fsaverage.BN_Atlas.32k_fs_LR.dlabel.nii';
glasser_L1 = ciftiopen(template_path1);
glasser_R1 = ciftiopen(template_path2);
cifti1 = ciftiopen(ciftipath);
path = 'E:\Phase dynamics and INT\Results\Visualization';

for iband = 1:2:8
    data = mean_phase_mean_stat2(:,iband);
    
    name = ['fp_band_',num2str(iband)];
    a_visualphasevar(data,cifti1,glasser_L1,glasser_R1,path,name);
end

