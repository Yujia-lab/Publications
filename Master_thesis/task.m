%% unzip file
clear
clc

path_rest = 'F:\HCP\FunRawCRSF\slow-5';
sub_file_rest = dir(path_rest);
sub_file_rest(1:2) = [];

for i = 1:328
    sub_id = sub_file_rest(i).name;
    sub_id = sub_id(1:6);
    SUB_ID(i) = str2num(sub_id);
end
SUB_ID = unique(SUB_ID);                  %select 82 subjects


path = 'F:\HCP_task\Zip_data\WM';
out_path = 'F:\HCP_task\WM';
sub_file = dir(path);
sub_file(1:2) = [];
sub_file(2:2:200) = [];
mkdir(out_path)


parfor isub = 1:length(sub_file)
    sub_path = fullfile(path,sub_file(isub).name);
    sub_id = sub_file(isub).name;
    sub_id = str2num(sub_id(1:6));
    if ismember(sub_id,SUB_ID) == 1
        unzip(sub_path,out_path)
    end
end

%% Move volume files
clc
clear
sub_path = 'F:\HCP_task\WM';
sub_file = dir(sub_path);
sub_file(1:2) = [];

out_path = 'F:\HCP_task\WM_vol';

for i = 1:length(sub_file)
    sub_id = sub_file(i).name;
    
    vol_path = fullfile(sub_path,sub_id,'MNINonLinear\Results\tfMRI_WM_LR');
    cd(vol_path);
    vol_name = dir('*LR.nii.gz');
    out_path2 = fullfile(out_path,[sub_id,'_LR']);
    mkdir(out_path2)
    movefile(fullfile(vol_path,vol_name.name),out_path2);           % move LR volume to outpath
    
    vol_path = fullfile(sub_path,sub_id,'MNINonLinear\Results\tfMRI_WM_RL');
    cd(vol_path);
    vol_name = dir('*RL.nii.gz');
    out_path2 = fullfile(out_path,[sub_id,'_RL']);
    mkdir(out_path2)
    movefile(fullfile(vol_path,vol_name.name),out_path2);           % move RL volume to outpath
end

%% Using dpabi to regress white matter and CSF
%% Move head motion and physiological files
clc
clear
sub_path = 'F:\HCP_task\WM';
sub_file = dir(sub_path);
sub_file(1:2) = [];

out_path1 = 'F:\HCP_task\WM_headmotion';
out_path2 = 'F:\HCP_task\WM_physio';
mkdir(out_path1)
mkdir(out_path2)

for i = 1:length(sub_file)
    sub_id = sub_file(i).name;
    
    hm_path = fullfile(sub_path,sub_id,'MNINonLinear\Results\tfMRI_WM_LR\Movement_Regressors_dt.txt');
    phy_path = fullfile(sub_path,sub_id,'MNINonLinear\Results\tfMRI_WM_LR\tfMRI_WM_LR_Physio_log.txt');
    copyfile(hm_path,out_path1);           % move headmotion parameters to outpath1
    copyfile(phy_path,out_path2);           % move physiological parameters to outpath2
    cd(out_path1)
    oldname = 'Movement_Regressors_dt.txt';
    newname = [sub_id,'_LR.txt'];
    eval(['!ren' 32 oldname 32 newname]);      % rename txt
    cd(out_path2)
    oldname = 'tfMRI_WM_LR_Physio_log.txt';
    newname = [sub_id,'_LR.txt'];
    eval(['!ren' 32 oldname 32 newname]);
    
    
    hm_path = fullfile(sub_path,sub_id,'MNINonLinear\Results\tfMRI_WM_RL\Movement_Regressors_dt.txt');
    phy_path = fullfile(sub_path,sub_id,'MNINonLinear\Results\tfMRI_WM_RL\tfMRI_WM_RL_Physio_log.txt');
    copyfile(hm_path,out_path1);           % move headmotion parameters to outpath1
    copyfile(phy_path,out_path2);           % move physiological parameters to outpath2
    cd(out_path1)
    oldname = 'Movement_Regressors_dt.txt';
    newname = [sub_id,'_RL.txt'];
    eval(['!ren' 32 oldname 32 newname]);      % rename txt
    cd(out_path2)
    oldname = 'tfMRI_WM_RL_Physio_log.txt';
    newname = [sub_id,'_RL.txt'];
    eval(['!ren' 32 oldname 32 newname]);
end

%% downsample physiological regressors
clc
clear
phy_path = 'F:\HCP_task\WM_physio';
sub_file = dir(phy_path);
sub_file(1:2) = [];
out_path = 'F:\HCP_task\WM_physio_resample';
mkdir(out_path);


for i = 1:length(sub_file)
    cd(phy_path)
    phy_file = importdata(sub_file(i).name);
    phy1 = resample(phy_file(:,2),10,2885);
    phy2 = resample(phy_file(:,3),10,2885);
    phy_resample = [phy1,phy2];
    cd(out_path)
    sub_id = sub_file(i).name;
    sub_id = sub_id(1:9);
    save(sub_id,'phy_resample');
end

%% regressing physiological and headmotion factors
clc
clear
vol_path = 'F:\HCP_task\WM_volC';
hm_path = 'F:\HCP_task\WM_headmotion';
phy_path = 'F:\HCP_task\WM_physio_resample';
out_path = 'F:\HCP_task\WM_volCR';

vol_file = dir(vol_path);
hm_file = dir(hm_path);
phy_file = dir(phy_path);
vol_file(1:2) = [];
hm_file(1:2) = [];
phy_file(1:2) = [];

for isubj = 1:length(vol_file)
    isubj
    hm = importdata(fullfile(hm_path,hm_file(isubj).name));
    load(fullfile(phy_path,phy_file(isubj).name));
    
    [volume, header] = y_Read(fullfile(vol_path, vol_file(isubj).name, 'CovRegressed_4DVolume.nii'));
    dim = size(volume);
    volume2D = reshape(volume, dim(1)*dim(2)*dim(3),dim(4));
    x1 = [hm, ones(dim(4),1)];
    x2 = [phy_resample, ones(dim(4),1)];
    out_volume = zeros(dim(1)*dim(2)*dim(3),dim(4));
    
    for ivoxel = 1:dim(1)*dim(2)*dim(3);
        if max(volume2D(ivoxel,:)) - min(volume2D(ivoxel,:)) ~= 0
            [b,bint,r,rint,stats] = regress(volume2D(ivoxel,:)',x1);
            [b,bint,r2,rint,stadts] = regress(r,x2);
            out_volume(ivoxel,:) = r2;
        end
    end
    
    out_volume = reshape(out_volume, dim(1),dim(2),dim(3),dim(4));
    mkdir(fullfile(out_path,vol_file(isubj).name));
    cd(fullfile(out_path,vol_file(isubj).name));
    y_Write(out_volume, header, 'CovRegressed_4DVolume.nii');
end

%% bands filter
clear
clc
bands = ccs_core_lfobands(405, 0.72);
for iband = 1:5
    bands{iband} = roundn(bands{iband},-4);
end

bands{5} = [0.6065,0.6844];
sub_path = 'F:\HCP_task\WM_volCR';
out_path = 'F:\HCP_task\WM_volCRF';
Fs = 1/0.72;

a_Nii_FIRfilter(sub_path,out_path,bands,Fs);

%% Extracting GS and ROI signals
clear
clc
sub_paths{1} = 'F:\HCP_task\WM_volCRF\slow-1';
sub_paths{2} = 'F:\HCP_task\WM_volCRF\slow-2';
sub_paths{3} = 'F:\HCP_task\WM_volCRF\slow-3';
sub_paths{4} = 'F:\HCP_task\WM_volCRF\slow-4';
sub_paths{5} = 'F:\HCP_task\WM_volCRF\slow-5';

out_paths{1} = 'E:\HCP\Results_task\slow-1';
out_paths{2} = 'E:\HCP\Results_task\slow-2';
out_paths{3} = 'E:\HCP\Results_task\slow-3';
out_paths{4} = 'E:\HCP\Results_task\slow-4';
out_paths{5} = 'E:\HCP\Results_task\slow-5';

mask_path = 'D:\program\matlab\Common\Mask\BN_Atlas_246_2mm.nii';

for i = 1:5
    [ROISignals,GS] = a_GSROI(sub_paths{i},mask_path);
    mkdir(out_paths{i})
    cd(out_paths{i});
    save ROISignals ROISignals
    save GS GS
end

%% Compute GS topography in each phase (fcmat是两个run合起来算的）
clear
clc

sub_paths{1} = 'E:\HCP\Results_task\slow-1';
sub_paths{2} = 'E:\HCP\Results_task\slow-2';
sub_paths{3} = 'E:\HCP\Results_task\slow-3';
sub_paths{4} = 'E:\HCP\Results_task\slow-4';
sub_paths{5} = 'E:\HCP\Results_task\slow-5';
% GStopo and FCmat
for i = 1:5
    i
    load(fullfile(sub_paths{i},'GS.mat'))
    load(fullfile(sub_paths{i},'ROISignals.mat'))
    GS_imag = hilbert(GS);
    GS_phase = angle(GS_imag);
    GS_phase = GS_phase * 180 / pi;
    GS = squeeze(GS);
    cd(sub_paths{i});
    save GS_phase GS_phase
    %sldwd_GStopo(GS,GS_phase,ROISignals,30,1,sub_paths{i},0)
    sldwd_fcmat_task(GS,GS_phase,ROISignals,30,1,sub_paths{i},0)
end

%% display strength and correlation
sub_paths{1} = 'E:\HCP\Results_task\slow-1';
sub_paths{2} = 'E:\HCP\Results_task\slow-2';
sub_paths{3} = 'E:\HCP\Results_task\slow-3';
sub_paths{4} = 'E:\HCP\Results_task\slow-4';
sub_paths{5} = 'E:\HCP\Results_task\slow-5';

for i = 1:5
    path = fullfile(sub_paths{i},'GS_topo_mean.mat');
    load(path)
    GS_topo_mean = [GS_topo_mean(:,end-14:end),GS_topo_mean(:,1:end-15)];
    CVC = corr(GS_topo_mean,GS_topo_mean,'type','pearson');
    cd(sub_paths{i})
    save CVC CVC
end

%% Phase difference of GStopo and GS phase
clear
clc

sub_paths{1} = 'E:\HCP\Results_task\slow-1';
sub_paths{2} = 'E:\HCP\Results_task\slow-2';
sub_paths{3} = 'E:\HCP\Results_task\slow-3';
sub_paths{4} = 'E:\HCP\Results_task\slow-4';
sub_paths{5} = 'E:\HCP\Results_task\slow-5';

for i = 1:5
    load(fullfile(sub_paths{i},'GS_topo_mean.mat'))
    GS_topo_mean = [GS_topo_mean(:,end-14:end),GS_topo_mean(:,1:end-15)];
    GS_topo_mean = GS_topo_mean';
    GS_topo_mean_hilb = hilbert(GS_topo_mean);
    strength_phase = angle(GS_topo_mean_hilb);
    strength_phase = strength_phase * 180 / pi;
    
    GS_topo_mean_all = mean(GS_topo_mean,2);
    GS_topo_mean_all_hilb = hilbert(GS_topo_mean_all);
    strength_phase_all = angle(GS_topo_mean_all_hilb);
    strength_phase_all = strength_phase_all * 180 / pi;
    
    phase_difference = bsxfun(@minus,strength_phase,strength_phase_all);
    CVC_phase = corr(phase_difference',phase_difference','type','pearson');
    cd(sub_paths{i});
    save intensity_phase strength_phase phase_difference CVC_phase
end

%% 网络效率
clear
clc

sub_paths{1} = 'E:\HCP\Results_task\slow-1';
sub_paths{2} = 'E:\HCP\Results_task\slow-2';
sub_paths{3} = 'E:\HCP\Results_task\slow-3';
sub_paths{4} = 'E:\HCP\Results_task\slow-4';
sub_paths{5} = 'E:\HCP\Results_task\slow-5';

phase_range = [-180:179];

for i = [3,4,5]
    i
    load(fullfile(sub_paths{i},'intensity_phase.mat'));
    
    for icl = 1:360
        x = CVC_phase(:,icl);
        x(icl) = [];
        mean_cl_z(icl) = mean(abs(a_fishertrans(x)));
    end
    mean_cl_r = gretna_inv_fishertrans(mean_cl_z);
    cd(sub_paths{i})
    
    mean_cl_p = a_r2p(mean_cl_r,246).*360;
    x = mean_cl_r;
    x(mean_cl_p > 0.05) = 1;
    min_r(i) = min(abs(x));
    save mean_cl_r mean_cl_r
    
    
    range1 = [35:75];
    range2 = [90:116];
    range3 = [215:265];
    range4 = [266:360];
    [min_value,pos1] = min(mean_cl_r(range1));
    TA_phase(i) = -180 + (pos1 + 34);
    [min_value,pos2] = min(mean_cl_r(range2));
    AP_phase(i) = -180 + (pos2 + 89);
    [min_value,pos3] = min(mean_cl_r(range3));
    PD_phase(i) = -180 + (pos3 + 214);
    [min_value,pos4] = min(mean_cl_r(range4));
    DT_phase(i) = -180 + (pos4 + 265);
    trough = [-180:TA_phase(i),DT_phase(i):179];
    peak = [AP_phase(i):PD_phase(i)];
    ascend = [TA_phase(i):AP_phase(i)];
    descend = [PD_phase(i):DT_phase(i)];
    trough_bin = ismember(phase_range,trough);
    peak_bin = ismember(phase_range,peak);
    ascend_bin = ismember(phase_range,ascend);
    descend_bin = ismember(phase_range,descend);
end


sub_fcmat_path = fullfile(sub_paths{i},'fcmat');
fcmat_file = dir(sub_fcmat_path);
fcmat_file(1:2) = [];

for isub = 1:length(fcmat_file)
    load(fullfile(sub_fcmat_path,fcmat_file(isub).name));
    fcmat = cat(3,fcmat(:,:,end-14:end),fcmat(:,:,1:end-15));
    fcmat = a_fishertrans(fcmat);   %fisher变换以做不同R值的运算
    trough_mat = mean(fcmat(:,:,trough_bin),3);    %平均波谷
    trough_mat = gretna_inv_fishertrans(trough_mat);     % fisher逆变换以计算网络效率
    nanpos = isnan(trough_mat);
    trough_mat(nanpos) = 0;
    K_value = prctile(reshape(trough_mat,[],1),75);
    trough_mat(trough_mat < K_value) = 0;     %小于0的R值置为0
    El_trough(:,isub) = efficiency_wei(trough_mat, 2);
    Eg_trough(isub) = efficiency_wei(trough_mat);
    
    ascend_mat = mean(fcmat(:,:,ascend_bin),3);    %平均波谷
    ascend_mat = gretna_inv_fishertrans(ascend_mat);     % fisher逆变换以计算网络效率
    nanpos = isnan(ascend_mat);
    ascend_mat(nanpos) = 0;
    K_value = prctile(reshape(ascend_mat,[],1),75);
    ascend_mat(ascend_mat < K_value) = 0;     %小于0的R值置为0
    El_ascend(:,isub) = efficiency_wei(ascend_mat, 2);
    Eg_ascend(isub) = efficiency_wei(ascend_mat);
    
    peak_mat = mean(fcmat(:,:,peak_bin),3);    %平均波谷
    peak_mat = gretna_inv_fishertrans(peak_mat);     % fisher逆变换以计算网络效率
    nanpos = isnan(peak_mat);
    peak_mat(nanpos) = 0;
    K_value = prctile(reshape(peak_mat,[],1),75);
    peak_mat(peak_mat < K_value) = 0;     %小于0的R值置为0
    El_peak(:,isub) = efficiency_wei(peak_mat, 2);
    Eg_peak(isub) = efficiency_wei(peak_mat);
    
    descend_mat = mean(fcmat(:,:,descend_bin),3);    %平均波谷
    descend_mat = gretna_inv_fishertrans(descend_mat);     % fisher逆变换以计算网络效率
    nanpos = isnan(descend_mat);
    descend_mat(nanpos) = 0;
    K_value = prctile(reshape(descend_mat,[],1),75);
    descend_mat(descend_mat < K_value) = 0;     %小于0的R值置为0
    El_descend(:,isub) = efficiency_wei(descend_mat, 2);
    Eg_descend(isub) = efficiency_wei(descend_mat);
end
cd(sub_paths{i})
save Efficiency Eg_trough Eg_ascend Eg_peak Eg_descend El_trough El_ascend El_peak El_descend
%% 四种状态网络效率的方差分析
clear
clc

sub_paths{1} = 'E:\HCP\Results_task\slow-1';
sub_paths{2} = 'E:\HCP\Results_task\slow-2';
sub_paths{3} = 'E:\HCP\Results_task\slow-3';
sub_paths{4} = 'E:\HCP\Results_task\slow-4';
sub_paths{5} = 'E:\HCP\Results_task\slow-5';

for i = 1:5

    load(fullfile(sub_paths{i},'Efficiency'));
    t = table(Eg_trough',Eg_ascend',Eg_peak',Eg_descend',...
        'VariableNames',{'y1','y2','y3','y4'});
    meas = table({'T','A','P','D'}','VariableNames',{'phase'}');
    rm = fitrm(t,'y1-y4~1','WithinDesign',meas);
    [tbl,A,C,D] = ranova(rm,'WithinModel','phase');
    F(i) = table2array(tbl(3,4));
    P(i) = table2array(tbl(3,6));
end
    eta_sq = a_etasq(F,3,243);

    
    for iROI = 1:246
        t = table(El_trough(iROI,:)',El_ascend(iROI,:)',El_peak(iROI,:)',El_descend(iROI,:)',...
            'VariableNames',{'y1','y2','y3','y4'});
        meas = table({'T','A','P','D'}','VariableNames',{'phase'}');
        rm = fitrm(t,'y1-y4~1','WithinDesign',meas);
        [tbl,A,C,D] = ranova(rm,'WithinModel','phase');
        P_anova(iROI,i) = table2array(tbl(3,6)).*246;
        F(iROI,i) = table2array(tbl(3,4));
    end
    F(P_anova(:,i) > 0.05,i) = 0;
    cd(sub_paths{i});
    %     a_vec2nii(F(:,i),'EL_anova.nii')
    %     a_vec2nii(El_trough(:,i),'EL_trough.nii')
    %     a_vec2nii(El_ascend(:,i),'El_ascend.nii')
    %     a_vec2nii(El_peak(:,i),'El_peak.nii')
    %     a_vec2nii(El_descend(:,i),'El_descend.nii')
    
    % 单T检验
    
    mean_El_ascend(:,i) = mean(El_ascend,2);
    mean_El_peak(:,i) = mean(El_peak,2);
    mean_El_descend(:,i) = mean(El_descend,2);
    mean_El_trough(:,i) = mean(El_trough,2);
    for iROI = 1:246
        [h,p,ci,stats] = ttest(El_ascend(iROI,:));
        El_ascend_T(iROI,i) = stats.tstat;
        [h,p,ci,stats] = ttest(El_peak(iROI,:));
        El_peak_T(iROI,i) = stats.tstat;
        [h,p,ci,stats] = ttest(El_descend(iROI,:));
        El_descend_T(iROI,i) = stats.tstat;
        [h,p,ci,stats] = ttest(El_trough(iROI,:));
        El_trough_T(iROI,i) = stats.tstat;
    end
    
     for i1 = 1:246
        for i2 = 1:246
            [h,P_descend(i1,i2),ci,stats] = ttest(El_descend(i1,:),El_descend(i2,:));
            T_descend(i1,i2) = stats.tstat;
        end
    end
    FWE_T_descend = a_multicorrect(T_descend,P_descend,1);
    FWE_T_descend(isnan(FWE_T_descend)==1) = 0;
    
    for i1 = 1:246
        for i2 = 1:246
            [h,P_ascend(i1,i2),ci,stats] = ttest(El_ascend(i1,:),El_ascend(i2,:));
            T_ascend(i1,i2) = stats.tstat;
        end
    end
    FWE_T_ascend = a_multicorrect(T_ascend,P_ascend,1);
    FWE_T_ascend(isnan(FWE_T_ascend)==1) = 0;
    
    for i1 = 1:246
        for i2 = 1:246
            [h,P_trough(i1,i2),ci,stats] = ttest(El_trough(i1,:),El_trough(i2,:));
            T_trough(i1,i2) = stats.tstat;
        end
    end
    FWE_T_trough = a_multicorrect(T_trough,P_trough,1);
    FWE_T_trough(isnan(FWE_T_trough)==1) = 0;
    
    for i1 = 1:246
        for i2 = 1:246
            [h,P_peak(i1,i2),ci,stats] = ttest(El_peak(i1,:),El_peak(i2,:));
            T_peak(i1,i2) = stats.tstat;
        end
    end
    FWE_T_peak = a_multicorrect(T_peak,P_peak,1);
    FWE_T_peak(isnan(FWE_T_peak)==1) = 0;
    
    cd(sub_paths{i})   
    save FWE_T FWE_T_peak FWE_T_trough FWE_T_ascend FWE_T_descend
    save El_T El_ascend_T El_peak_T El_descend_T El_trough_T F
    save mean_El mean_El_ascend mean_El_peak mean_El_descend mean_El_trough F
end


for i1 = 1:246
    for i2 = 1:246
        [h,p(i1,i2),ci,stats] = ttest(El_descend(i1,:),El_descend(i2,:));
        T(i1,i2) = stats.tstat;
    end
end


%% sorting correct behavioral data
clc
clear
% import txt data
sub_path = 'F:\HCP_task\WM';
sub_file = dir(sub_path);
sub_file(1:2) = [];

for i = 1:length(sub_file)
    i
    TAB_path1 = fullfile(sub_path,sub_file(i).name,'MNINonLinear\Results\tfMRI_WM_LR');
    % sort RT
    cd(TAB_path1);
    file1 = dir('WM*TAB.txt');
    TAB1 = importdata(file1.name);
    TAB1 = TAB1.textdata;
    stim_onset = a_cell2mat(TAB1(:,87));
    stim_RTtime = a_cell2mat(TAB1(:,89));
    stim_ACC = a_cell2mat(TAB1(:,90));
    stim_RT = a_cell2mat(TAB1(:,91));
    stim_con = a_cell2mat(TAB1(:,74));
    stim_info{2*i-1} = [stim_con,stim_onset,stim_RTtime,stim_RT,stim_ACC];
    % sort onset
    cd(fullfile(TAB_path1,'EVs'))
    file_0bkcor = dir('0bk_cor.txt');
    file_2bkcor = dir('2bk_cor.txt');
    bk0cor_trial = importdata(file_0bkcor.name);
    bk2cor_trial = importdata(file_2bkcor.name);
    bk0cor_onsets{2*i-1} = bk0cor_trial;
    bk2cor_onsets{2*i-1} = bk2cor_trial;
    
    % sort RT
    TAB_path2 = fullfile(sub_path,sub_file(i).name,'MNINonLinear\Results\tfMRI_WM_RL');
    cd(TAB_path2);
    file2 = dir('WM*TAB.txt');
    TAB2 = importdata(file2.name);
    TAB2 = TAB2.textdata;
    stim_onset = a_cell2mat(TAB2(:,87));
    stim_RTtime = a_cell2mat(TAB2(:,89));
    stim_ACC = a_cell2mat(TAB2(:,90));
    stim_RT = a_cell2mat(TAB2(:,91));
    stim_con = a_cell2mat(TAB2(:,74));
    stim_info{2*i} = [stim_con,stim_onset,stim_RTtime,stim_RT,stim_ACC];
    % sort onset
    cd(fullfile(TAB_path2,'EVs'))
    file_0bkcor = dir('0bk_cor.txt');
    file_2bkcor = dir('2bk_cor.txt');
    bk0cor_trial = importdata(file_0bkcor.name);
    bk2cor_trial = importdata(file_2bkcor.name);
    bk0cor_onsets{2*i} = bk0cor_trial;
    bk2cor_onsets{2*i} = bk2cor_trial;
end

cd('E:\HCP\Results_task')
save behav_data_cor stim_info bk0cor_onsets bk2cor_onsets

%% sorting err behavioral data
clc
clear
% import txt data
sub_path = 'F:\HCP_task\WM';
sub_file = dir(sub_path);
sub_file(1:2) = [];

for i = 1:length(sub_file)
    i
    TAB_path1 = fullfile(sub_path,sub_file(i).name,'MNINonLinear\Results\tfMRI_WM_LR');
    % sort RT
    cd(TAB_path1);
    file1 = dir('WM*TAB.txt');
    TAB1 = importdata(file1.name);
    TAB1 = TAB1.textdata;
    stim_onset = a_cell2mat(TAB1(:,87));
    stim_RTtime = a_cell2mat(TAB1(:,89));
    stim_ACC = a_cell2mat(TAB1(:,90));
    stim_RT = a_cell2mat(TAB1(:,91));
    stim_con = a_cell2mat(TAB1(:,74));
    stim_info{2*i-1} = [stim_con,stim_onset,stim_RTtime,stim_RT,stim_ACC];
    % sort onset
    cd(fullfile(TAB_path1,'EVs'))
    file_0bkerr = dir('0bk_err.txt');
    file_2bkerr = dir('2bk_err.txt');
    bk0err_trial = importdata(file_0bkerr.name);
    bk2err_trial = importdata(file_2bkerr.name);
    bk0err_onsets{2*i-1} = bk0err_trial;
    bk2err_onsets{2*i-1} = bk2err_trial;
    
    % sort RT
    TAB_path2 = fullfile(sub_path,sub_file(i).name,'MNINonLinear\Results\tfMRI_WM_RL');
    cd(TAB_path2);
    file2 = dir('WM*TAB.txt');
    TAB2 = importdata(file2.name);
    TAB2 = TAB2.textdata;
    stim_onset = a_cell2mat(TAB2(:,87));
    stim_RTtime = a_cell2mat(TAB2(:,89));
    stim_ACC = a_cell2mat(TAB2(:,90));
    stim_RT = a_cell2mat(TAB2(:,91));
    stim_con = a_cell2mat(TAB2(:,74));
    stim_info{2*i} = [stim_con,stim_onset,stim_RTtime,stim_RT,stim_ACC];
    % sort onset
    cd(fullfile(TAB_path2,'EVs'))
    file_0bkerr = dir('0bk_err.txt');
    file_2bkerr = dir('2bk_err.txt');
    bk0err_trial = importdata(file_0bkerr.name);
    bk2err_trial = importdata(file_2bkerr.name);
    bk0err_onsets{2*i} = bk0err_trial;
    bk2err_onsets{2*i} = bk2err_trial;
end

cd('E:\HCP\Results_task')
save behav_data_err stim_info bk0err_onsets bk2err_onsets



%% 0-back reaction time and GS phase
clear
clc
sub_paths{1} = 'E:\HCP\Results_task\slow-1';
sub_paths{2} = 'E:\HCP\Results_task\slow-2';
sub_paths{3} = 'E:\HCP\Results_task\slow-3';
sub_paths{4} = 'E:\HCP\Results_task\slow-4';
sub_paths{5} = 'E:\HCP\Results_task\slow-5';

load('E:\HCP\Results_task\behav_data.mat')
load('E:\HCP\Results_task\behav_data_err.mat')


time = [0:0.72:404*0.72];
lag = [0:0.5:10];

for i = 1:5
    load(fullfile(sub_paths{i},'GS_phase.mat'));
    RT_all = [];
    phase_interp_all = [];
    RT_norm_all = [];
    for isub = 1:length(stim_info)
        
        %         isub = 1 + (isub-1)*2;              %LR的被试
        sub_phase = GS_phase(:,isub);
        sub_behav = stim_info{isub};
        
        %分离0-back的RT
        
        bk0_RT = sub_behav(sub_behav(:,1)==0,:);
        
        %提取正确0-back的onsets
        bk0_onsets = bk0cor_onsets{isub};
        
        % 插值正确0-back的GS相位
        for itrial = 1:size(bk0_onsets,1)
            for ilag = 1:length(lag)
                phase_interp = interp1(time',sub_phase,bk0_onsets(itrial,1)+lag(ilag),'linear');
                if phase_interp < -180
                    phase_interp = phase_interp + 360
                elseif phase_interp > 180
                    phase_interp = phase_interp - 360
                end
                bk0_phase_interp_all_cor{isub}(itrial,ilag) = phase_interp;         %所有被试的相位
            end
        end
        
        %提取错误0-back的onsets
        bk0_onsets = bk0err_onsets{isub};
        % 插值错误0-back的GS相位
        for itrial = 1:size(bk0_onsets,1)
            for ilag = 1:length(lag)
                phase_interp = interp1(time',sub_phase,bk0_onsets(itrial,1)+lag(ilag),'linear');
                if phase_interp < -180
                    phase_interp = phase_interp + 360
                elseif phase_interp > 180
                    phase_interp = phase_interp - 360
                end
                bk0_phase_interp_all_err{isub}(itrial,ilag) = phase_interp;         %所有被试的相位
            end
        end
        
        
        bk0_RT_all{isub} = [bk0_RT(:,4),bk0_RT(:,5)];         %所有被试的反应时和正确率
    end
    
    cd(sub_paths{i})
    save bk0_info bk0_phase_interp_all_cor bk0_phase_interp_all_err bk0_RT_all
end

%% 2-back reaction time and GS phase
clear
clc
sub_paths{1} = 'E:\HCP\Results_task\slow-1';
sub_paths{2} = 'E:\HCP\Results_task\slow-2';
sub_paths{3} = 'E:\HCP\Results_task\slow-3';
sub_paths{4} = 'E:\HCP\Results_task\slow-4';
sub_paths{5} = 'E:\HCP\Results_task\slow-5';

load('E:\HCP\Results_task\behav_data.mat')
load('E:\HCP\Results_task\behav_data_err.mat')

time = [0:0.72:404*0.72];
lag = [0:0.5:10];

for i = 1:5
    load(fullfile(sub_paths{i},'GS_phase.mat'));
    RT_all = [];
    phase_interp_all = [];
    RT_norm_all = [];
    for isub = 1:length(stim_info)
        
        %         isub = 1 + (isub-1)*2;              %LR的被试
        sub_phase = GS_phase(:,isub);
        sub_behav = stim_info{isub};
        
        %分离2-back的RT
        bk2_RT = sub_behav(sub_behav(:,1)==2,:);
        
        %提取正确2-back的onsets
        bk2_onsets = bk2cor_onsets{isub};
        
        % 插值正确2-back的GS相位
        for itrial = 1:size(bk2_onsets,1)
            for ilag = 1:length(lag)
                phase_interp = interp1(time',sub_phase,bk2_onsets(itrial,1)+lag(ilag),'linear');
                if phase_interp < -180
                    phase_interp = phase_interp + 360
                elseif phase_interp > 180
                    phase_interp = phase_interp - 360
                end
                bk2_phase_interp_all_cor{isub}(itrial,ilag) = phase_interp;         %所有被试的相位
            end
        end
        
        %提取错误2-back的onsets
        bk2_onsets = bk2err_onsets{isub};
        % 插值错误2-back的GS相位
        for itrial = 1:size(bk2_onsets,1)
            for ilag = 1:length(lag)
                phase_interp = interp1(time',sub_phase,bk2_onsets(itrial,1)+lag(ilag),'linear');
                if phase_interp < -180
                    phase_interp = phase_interp + 360
                elseif phase_interp > 180
                    phase_interp = phase_interp - 360
                end
                bk2_phase_interp_all_err{isub}(itrial,ilag) = phase_interp;         %所有被试的相位
            end
        end
        
        bk2_RT_all{isub} = [bk2_RT(:,4),bk2_RT(:,5)];         %所有被试的反应时
        
    end
    
    cd(sub_paths{i})
    save bk2_info bk2_phase_interp_all_cor bk2_phase_interp_all_err bk2_RT_all
end

%% 0-back
clear
clc
sub_paths{1} = 'E:\HCP\Results_task\slow-1';
sub_paths{2} = 'E:\HCP\Results_task\slow-2';
sub_paths{3} = 'E:\HCP\Results_task\slow-3';
sub_paths{4} = 'E:\HCP\Results_task\slow-4';
sub_paths{5} = 'E:\HCP\Results_task\slow-5';



for i = 1:5
    load(fullfile(sub_paths{i},'bk0_info.mat'));
    
    load(fullfile(sub_paths{i},'intensity_phase'));
    mean_CVC = mean(CVC_phase);
    mean_CVC_smooth = smooth(mean_CVC);
    mean_CVC_smooth_abs = abs(mean_CVC_smooth);
    range1 = [35:75];
    range2 = [76:116];
    range3 = [215:265];
    range4 = [266:360];
    [min_value,pos1] = min(mean_CVC_smooth_abs(range1));
    TA_phase(i) = -180 + (pos1 + 34)+ 15;
    [min_value,pos2] = min(mean_CVC_smooth_abs(range2));
    AP_phase(i) = -180 + (pos2 + 75)+ 15;
    [min_value,pos3] = min(mean_CVC_smooth_abs(range3));
    PD_phase(i) = -180 + (pos3 + 214)+ 15;
    [min_value,pos4] = min(mean_CVC_smooth_abs(range4));
    DT_phase(i) = -180 + (pos4 + 265)+ 15;
    
    
    for isub = 1:length(bk0_RT_all)/2
        
        ISUB = isub;
        isub = 1 + (isub-1)*2;
        
        
        bk0_phase_LR = bk0_phase_interp_all_cor{isub};
        bk0_phase_RL = bk0_phase_interp_all_cor{isub+1};
        bk0_phase = [bk0_phase_LR;bk0_phase_RL];
        
        bk0_RT_LR = bk0_RT_all{isub};
        bk0_RT_RL = bk0_RT_all{isub+1};
        bk0_RT = [bk0_RT_LR;bk0_RT_RL];
        bk0_RT(bk0_RT(:,2) == 0,:) = []; %删除错误试次
        bk0_RT(:,2) = [];
        
        
        amp = bk0_RT; %使用原始反应时，不使用就删掉
        for ilag = 1:size(bk0_phase,2)
            phase = bk0_phase(:,ilag);
            
            trough = find((phase >= -180 & phase <= TA_phase(i)) | (phase >= DT_phase(i) & phase <= 180));
            peak = find(phase >= AP_phase(i) & phase <= PD_phase(i));
            ascend = find(phase >= TA_phase(i) & phase <= AP_phase(i));
            descend = find(phase >= PD_phase(i) & phase <= DT_phase(i));
            
            %                 trough = find((phase >= -180 & phase <= -135) | (phase >= 135 & phase <= 180));
            %                 peak = find(phase >= -45 & phase <= 45);
            %                 ascend = find(phase >= -135 & phase <= -45);
            %                 descend = find(phase >= 45 & phase <= 135);
            
            amp_trough(ilag) = mean(amp(trough));
            amp_peak(ilag) = mean(amp(peak));
            amp_ascend(ilag) = mean(amp(ascend));
            amp_descend(ilag) =mean(amp(descend));
        end
        amp_trough(isnan(amp_trough) == 1) = mean(amp_trough(isnan(amp_trough) == 0));
        amp_peak(isnan(amp_peak) == 1) = mean(amp_peak(isnan(amp_peak) == 0));
        amp_ascend(isnan(amp_ascend) == 1) = mean(amp_ascend(isnan(amp_ascend) == 0));
        amp_descend(isnan(amp_descend) == 1) = mean(amp_descend(isnan(amp_descend) == 0));
        T0(:,ISUB) = amp_trough;
        P0(:,ISUB) = amp_peak;
        A0(:,ISUB) = amp_ascend;
        D0(:,ISUB) = amp_descend;
    end
    for ilag = 1:size(bk0_phase,2)
        Trough = squeeze(T0(ilag,i,:));
        Peak = squeeze(P0(ilag,i,:));
        Ascend = squeeze(A0(ilag,i,:));
        Descend = squeeze(D0(ilag,i,:));
        [p{ilag,i},tbl,stats{ilag,i}] = anova1([Peak,Trough,Ascend,Descend],[],'off');
    end
    cd(sub_paths{i});
    save RT_phase0 T0 P0 A0 D0
end



%% 2-back
clear
clc
sub_paths{1} = 'E:\HCP\Results_task\slow-1';
sub_paths{2} = 'E:\HCP\Results_task\slow-2';
sub_paths{3} = 'E:\HCP\Results_task\slow-3';
sub_paths{4} = 'E:\HCP\Results_task\slow-4';
sub_paths{5} = 'E:\HCP\Results_task\slow-5';



for i = 1:5
    load(fullfile(sub_paths{i},'bk2_info.mat'));
    
    load(fullfile(sub_paths{i},'intensity_phase'));
    mean_CVC = mean(CVC_phase);
    mean_CVC_smooth = smooth(mean_CVC);
    mean_CVC_smooth_abs = abs(mean_CVC_smooth);
    range1 = [35:75];
    range2 = [76:116];
    range3 = [215:265];
    range4 = [266:360];
    [min_value,pos1] = min(mean_CVC_smooth_abs(range1));
    TA_phase(i) = -180 + (pos1 + 34)+ 15;
    [min_value,pos2] = min(mean_CVC_smooth_abs(range2));
    AP_phase(i) = -180 + (pos2 + 75)+ 15;
    [min_value,pos3] = min(mean_CVC_smooth_abs(range3));
    PD_phase(i) = -180 + (pos3 + 214)+ 15;
    [min_value,pos4] = min(mean_CVC_smooth_abs(range4));
    DT_phase(i) = -180 + (pos4 + 265)+ 15;
    
    
    for isub = 1:length(bk2_RT_all)/2
        ISUB = isub;
        isub = 1 + (isub-1)*2;
        bk2_phase_LR = bk2_phase_interp_all_cor{isub};
        bk2_phase_RL = bk2_phase_interp_all_cor{isub+1};
        bk2_phase = [bk2_phase_LR;bk2_phase_RL];
        
        bk2_RT_LR = bk2_RT_all{isub};
        bk2_RT_RL = bk2_RT_all{isub+1};
        bk2_RT = [bk2_RT_LR;bk2_RT_RL];
        bk2_RT(bk2_RT(:,2) == 0,:) = []; %删除错误试次
        bk2_RT(:,2) = [];
        
        
        amp = bk2_RT; %使用原始反应时，不使用就删掉
        for ilag = 1:size(bk2_phase,2)
            phase = bk2_phase(:,ilag);
            
            trough = find((phase >= -180 & phase <= TA_phase(i)) | (phase >= DT_phase(i) & phase <= 180));
            peak = find(phase >= AP_phase(i) & phase <= PD_phase(i));
            ascend = find(phase >= TA_phase(i) & phase <= AP_phase(i));
            descend = find(phase >= PD_phase(i) & phase <= DT_phase(i));
            
            %                 trough = find((phase >= -180 & phase <= -135) | (phase >= 135 & phase <= 180));
            %                 peak = find(phase >= -45 & phase <= 45);
            %                 ascend = find(phase >= -135 & phase <= -45);
            %                 descend = find(phase >= 45 & phase <= 135);
            
            amp_trough(ilag) = mean(amp(trough));
            amp_peak(ilag) = mean(amp(peak));
            amp_ascend(ilag) = mean(amp(ascend));
            amp_descend(ilag) =mean(amp(descend));
        end
        amp_trough(isnan(amp_trough) == 1) = mean(amp_trough(isnan(amp_trough) == 0));
        amp_peak(isnan(amp_peak) == 1) = mean(amp_peak(isnan(amp_peak) == 0));
        amp_ascend(isnan(amp_ascend) == 1) = mean(amp_ascend(isnan(amp_ascend) == 0));
        amp_descend(isnan(amp_descend) == 1) = mean(amp_descend(isnan(amp_descend) == 0));
        T2(:,ISUB) = amp_trough;
        P2(:,ISUB) = amp_peak;
        A2(:,ISUB) = amp_ascend;
        D2(:,ISUB) = amp_descend;
    end
    
    for ilag = 1:size(bk2_phase,2)
        Trough = squeeze(T2(ilag,i,:));
        Peak = squeeze(P2(ilag,i,:));
        Ascend = squeeze(A2(ilag,i,:));
        Descend = squeeze(D2(ilag,i,:));
        [p{ilag,i},tbl,stats{ilag,i}] = anova1([Peak,Trough,Ascend,Descend],[],'off');
    end
    cd(sub_paths{i});
    save RT_phase2 T2 P2 A2 D2
    
end


%% 2（任务）X4（相位）方差分析
clear
clc
sub_paths{1} = 'E:\HCP\Results_task\slow-1';
sub_paths{2} = 'E:\HCP\Results_task\slow-2';
sub_paths{3} = 'E:\HCP\Results_task\slow-3';
sub_paths{4} = 'E:\HCP\Results_task\slow-4';
sub_paths{5} = 'E:\HCP\Results_task\slow-5';


for i = 1:5
    load(fullfile(sub_paths{i},'RT_phase0.mat'));
    load(fullfile(sub_paths{i},'RT_phase2.mat'));
    
    ilag = 1;
    t = table(T0(ilag,:)',P0(ilag,:)',A0(ilag,:)',D0(ilag,:)',...
        T2(ilag,:)',P2(ilag,:)',A2(ilag,:)',D2(ilag,:)',...
        'VariableNames',...
        {'y1','y2','y3','y4','y5','y6','y7','y8'});
    meas = table({'T','P','A','D','T','P','A','D'}',...
        {'0-back','0-back','0-back','0-back','2-back','2-back','2-back','2-back'}',...
        'VariableNames',{'phase','task'}');
    rm = fitrm(t,'y1-y8~1','WithinDesign',meas);
    tbl_ep = epsilon(rm);
    [tbl{i},A,C,D] = ranova(rm,'WithinModel',...
        'phase+task+phase*task');
    
end

F = 1.108915303590583;
dfa = 3;
dfe = 243;
ETA = a_etasq(F,dfa,dfe)


% 计算效应量
for i = 1:8
    SSA = table2array(tbl((i-1)*2+1,1));
    SSE = table2array(tbl(i*2,1));
    etasq = (SSA/(SSA+SSE))^2;
    cohenf(i) = sqrt(etasq/(1-etasq));
end

p_phase = p_phase .* 21;
p_inter = p_inter .* 21;
p_task = p_task .* 21;
cd('E:\HCP\Results_task\ANOVA');
save ANOVA_p p_phase p_task p_inter




% 0-back的简单主效应
for i = 1:5
    load(fullfile(sub_paths{i},'RT_phase0.mat'));
    MT0 = mean(T0,2);
    MA0 = mean(A0,2);
    
    
    MP0 = mean(P0,2);
    MD0 = mean(D0,2);
    
    
    MRT0 = [MT0,MA0,MP0,MD0];
    cd(sub_paths{i})
    save MRT0 MRT0
end

% simple effects
for ilag = 1:21
    [h,p,ci,stats] = ttest(P0(ilag,:),T0(ilag,:));
    T_P0(1,ilag) = p*60;
    [h,p,ci,stats] = ttest(P0(ilag,:),A0(ilag,:));
    T_P0(2,ilag) = p*60;
    [h,p,ci,stats] = ttest(P0(ilag,:),D0(ilag,:));
    T_P0(3,ilag) = p*60;
    [h,p,ci,stats] = ttest(T0(ilag,:),A0(ilag,:));
    T_P0(4,ilag) = p*60;
    [h,p,ci,stats] = ttest(T0(ilag,:),D0(ilag,:));
    T_P0(5,ilag) = p*60;
    [h,p,ci,stats] = ttest(A0(ilag,:),D0(ilag,:));
    T_P0(6,ilag) = p*60;
end

Cohend = 8.08/sqrt(81);

% 2-back的简单主效应
for i = 1:5
    load(fullfile(sub_paths{i},'RT_phase2.mat'));
    MT2 = mean(T2,2);
    MA2 = mean(A2,2);
    MP2 = mean(P2,2);
    MD2 = mean(D2,2);
    MRT2 = [MT2,MA2,MP2,MD2];
    cd(sub_paths{i})
    save MRT2 MRT2
end

% simple effects
for ilag = 1:21
    [h,p,ci,stats] = ttest(T2(i,:),A2(i,:));
    T_P2(1,ilag) = p;
    [h,p,ci,stats] = ttest(T2(i,:),P2(i,:));
    T_P2(2,ilag) = p;
    [h,p,ci,stats] = ttest(T2(i,:),D2(i,:));
    T_P2(3,ilag) = p;
    [h,p,ci,stats] = ttest(A2(i,:),P2(i,:));
    T_P2(4,ilag) = p;
    [h,p,ci,stats] = ttest(A2(i,:),D2(i,:));
    T_P2(5,ilag) = p;
    [h,p,ci,stats] = ttest(P2(i,:),D2(i,:));
    T_P2(6,ilag) = p;
end

%% 分析0-back正确率
clear
clc
sub_paths{1} = 'E:\HCP\Results_task\slow-1';
sub_paths{2} = 'E:\HCP\Results_task\slow-2';
sub_paths{3} = 'E:\HCP\Results_task\slow-3';
sub_paths{4} = 'E:\HCP\Results_task\slow-4';
sub_paths{5} = 'E:\HCP\Results_task\slow-5';

for i = 1:5
    load(fullfile(sub_paths{i},'bk0_info.mat'));
    load(fullfile(sub_paths{i},'intensity_phase'));
    mean_CVC = mean(CVC_phase);
    mean_CVC_smooth = smooth(mean_CVC);
    mean_CVC_smooth_abs = abs(mean_CVC_smooth);
    range1 = [35:75];
    range2 = [76:116];
    range3 = [215:265];
    range4 = [266:360];
    [min_value,pos1] = min(mean_CVC_smooth_abs(range1));
    TA_phase(i) = -180 + (pos1 + 34)+ 15;
    [min_value,pos2] = min(mean_CVC_smooth_abs(range2));
    AP_phase(i) = -180 + (pos2 + 75)+ 15;
    [min_value,pos3] = min(mean_CVC_smooth_abs(range3));
    PD_phase(i) = -180 + (pos3 + 214)+ 15;
    [min_value,pos4] = min(mean_CVC_smooth_abs(range4));
    DT_phase(i) = -180 + (pos4 + 265)+ 15;
    
    
    for isub = 1:length(bk0_RT_all)/2
        
        ISUB = isub;
        isub = 1 + (isub-1)*2;
        bk0_phase_LR = bk0_phase_interp_all_err{isub};
        bk0_phase_RL = bk0_phase_interp_all_err{isub+1};
        bk0_phase = [bk0_phase_LR;bk0_phase_RL];
        
        
        for ilag = 1:size(bk0_phase,2)
            phase = bk0_phase(:,ilag);
            
            trough = find((phase >= -180 & phase <= TA_phase(i)) | (phase >= DT_phase(i) & phase <= 180));
            peak = find(phase >= AP_phase(i) & phase <= PD_phase(i));
            ascend = find(phase >= TA_phase(i) & phase <= AP_phase(i));
            descend = find(phase >= PD_phase(i) & phase <= DT_phase(i));
            
            T_freq(ilag) = length(trough);
            P_freq(ilag) = length(peak);
            A_freq(ilag) = length(ascend);
            D_freq(ilag) = length(descend);
            
        end
        T0(:,ISUB) = T_freq;
        P0(:,ISUB) = P_freq;
        A0(:,ISUB) = A_freq;
        D0(:,ISUB) = D_freq;
    end
    A0 = (80 - A0) ./ 80;
    D0 = (80 - D0) ./ 80;
    P0 = (80 - P0) ./ 80;
    T0 = (80 - T0) ./ 80;
    cd(sub_paths{i});
    save RT_phase_ACC0 T0 P0 A0 D0
end

%% 分析2-back正确率
clear
clc
sub_paths{1} = 'E:\HCP\Results_task\slow-1';
sub_paths{2} = 'E:\HCP\Results_task\slow-2';
sub_paths{3} = 'E:\HCP\Results_task\slow-3';
sub_paths{4} = 'E:\HCP\Results_task\slow-4';
sub_paths{5} = 'E:\HCP\Results_task\slow-5';

for i = 1:5
    load(fullfile(sub_paths{i},'bk2_info.mat'));
    
    load(fullfile(sub_paths{i},'intensity_phase'));
    mean_CVC = mean(CVC_phase);
    mean_CVC_smooth = smooth(mean_CVC);
    mean_CVC_smooth_abs = abs(mean_CVC_smooth);
    range1 = [35:75];
    range2 = [76:116];
    range3 = [215:265];
    range4 = [266:360];
    [min_value,pos1] = min(mean_CVC_smooth_abs(range1));
    TA_phase(i) = -180 + (pos1 + 34)+ 15;
    [min_value,pos2] = min(mean_CVC_smooth_abs(range2));
    AP_phase(i) = -180 + (pos2 + 75)+ 15;
    [min_value,pos3] = min(mean_CVC_smooth_abs(range3));
    PD_phase(i) = -180 + (pos3 + 214)+ 15;
    [min_value,pos4] = min(mean_CVC_smooth_abs(range4));
    DT_phase(i) = -180 + (pos4 + 265)+ 15;
    
    
    for isub = 1:length(bk2_RT_all)/2
        
        ISUB = isub;
        isub = 1 + (isub-1)*2;
        bk2_phase_LR = bk2_phase_interp_all_err{isub};
        bk2_phase_RL = bk2_phase_interp_all_err{isub+1};
        bk2_phase = [bk2_phase_LR;bk2_phase_RL];
        
        
        for ilag = 1:size(bk2_phase,2)
            phase = bk2_phase(:,ilag);
            
            trough = find((phase >= -180 & phase <= TA_phase(i)) | (phase >= DT_phase(i) & phase <= 180));
            peak = find(phase >= AP_phase(i) & phase <= PD_phase(i));
            ascend = find(phase >= TA_phase(i) & phase <= AP_phase(i));
            descend = find(phase >= PD_phase(i) & phase <= DT_phase(i));
            
            T_freq(ilag) = length(trough);
            P_freq(ilag) = length(peak);
            A_freq(ilag) = length(ascend);
            D_freq(ilag) = length(descend);
        end
        
        T2(:,ISUB) = T_freq;
        P2(:,ISUB) = P_freq;
        A2(:,ISUB) = A_freq;
        D2(:,ISUB) = D_freq;
    end
    A2 = (80 - A2) ./ 80;
    D2 = (80 - D2) ./ 80;
    P2 = (80 - P2) ./ 80;
    T2 = (80 - T2) ./ 80;
    cd(sub_paths{i});
    save RT_phase_ACC2 T2 P2 A2 D2
end

%% 2（任务）X4（相位）方差分析 正确率
clear
clc
sub_paths{1} = 'E:\HCP\Results_task\slow-1';
sub_paths{2} = 'E:\HCP\Results_task\slow-2';
sub_paths{3} = 'E:\HCP\Results_task\slow-3';
sub_paths{4} = 'E:\HCP\Results_task\slow-4';
sub_paths{5} = 'E:\HCP\Results_task\slow-5';



for i = 1:5
    load(fullfile(sub_paths{i},'RT_phase_ACC0.mat'));
    load(fullfile(sub_paths{i},'RT_phase_ACC2.mat'));
    
    ilag = 1;
    t = table(T0(ilag,:)',P0(ilag,:)',A0(ilag,:)',D0(ilag,:)',...
        T2(ilag,:)',P2(ilag,:)',A2(ilag,:)',D2(ilag,:)',...
        'VariableNames',...
        {'y1','y2','y3','y4','y5','y6','y7','y8'});
    meas = table({'T','P','A','D','T','P','A','D'}',...
        {'0-back','0-back','0-back','0-back','2-back','2-back','2-back','2-back'}',...
        'VariableNames',{'phase','task'}');
    rm = fitrm(t,'y1-y8~1','WithinDesign',meas);
    tbl_ep = epsilon(rm);
    [tbl{i},A,C,D] = ranova(rm,'WithinModel',...
        'phase+task+phase*task');
    
end

F = 37.532669317804750;
dfa = 1;
dfe = 81;
ETA = a_etasq(F,dfa,dfe)


for i = 1:8
    SSA = table2array(tbl((i-1)*2+1,1));
    SSE = table2array(tbl(i*2,1));
    etasq = (SSA/(SSA+SSE))^2;
    cohenf(i) = sqrt(etasq/(1-etasq));
end
% bk0简单主效应

% 0-back的简单主效应
for i = 1:5
    load(fullfile(sub_paths{i},'RT_phase_ACC0.mat'));
    MT0 = mean(T0,2);
    MA0 = mean(A0,2);
    MP0 = mean(P0,2);
    MD0 = mean(D0,2);
    MRT0 = [MT0,MA0,MP0,MD0];
    cd(sub_paths{i})
    save MRT0_ACC MRT0
end

for ilag = 1:21
    [h,p,ci,stats] = ttest(P0(ilag,:),T0(ilag,:));
    T_P0(1,ilag) = p*60;
    [h,p,ci,stats] = ttest(P0(ilag,:),A0(ilag,:));
    T_P0(2,ilag) = p*60;
    [h,p,ci,stats] = ttest(P0(ilag,:),D0(ilag,:));
    T_P0(3,ilag) = p*60;
    [h,p,ci,stats] = ttest(T0(ilag,:),A0(ilag,:));
    T_P0(4,ilag) = p*60;
    [h,p,ci,stats] = ttest(T0(ilag,:),D0(ilag,:));
    T_P0(5,ilag) = p*60;
    [h,p,ci,stats] = ttest(A0(ilag,:),D0(ilag,:));
    T_P0(6,ilag) = p*60;
end


% bk2简单主效应
for i = 1:5
    load(fullfile(sub_paths{i},'RT_phase_ACC2.mat'));
    MT2 = mean(T2,2);
    MA2 = mean(A2,2);
    MP2 = mean(P2,2);
    MD2 = mean(D2,2);
    MRT2 = [MT2,MA2,MP2,MD2];
    cd(sub_paths{i})
    save MRT2_ACC MRT2
end
for ilag = 1:21
    [h,p,ci,stats] = ttest(T2(ilag,:),A2(ilag,:));
    T_P2(1,ilag) = p*60;
    [h,p,ci,stats] = ttest(T2(ilag,:),P2(ilag,:));
    T_P2(2,ilag) = p*60;
    [h,p,ci,stats] = ttest(T2(ilag,:),D2(ilag,:));
    T_P2(3,ilag) = p*60;
    [h,p,ci,stats] = ttest(A2(ilag,:),P2(ilag,:));
    T_P2(4,ilag) = p*60;
    [h,p,ci,stats] = ttest(A2(ilag,:),D2(ilag,:));
    T_P2(5,ilag) = p*60;
    [h,p,ci,stats] = ttest(P2(ilag,:),D2(ilag,:));
    T_P2(6,ilag) = p*60;
end


p_phase = p_phase .* 21;
p_inter = p_inter .* 21;
p_task = p_task .* 21;
cd('E:\HCP\Results_task\ANOVA');
save ANOVA_p_ACC p_phase p_task p_inter

%% global efficiency
clear
clc
sub_paths{1} = 'E:\HCP\Results_task\slow-1';
sub_paths{2} = 'E:\HCP\Results_task\slow-2';
sub_paths{3} = 'E:\HCP\Results_task\slow-3';
sub_paths{4} = 'E:\HCP\Results_task\slow-4';
sub_paths{5} = 'E:\HCP\Results_task\slow-5';

phase_range = [-165:180,-179:-166];

for i = 1:5
    i
    load(fullfile(sub_paths{i},'bk2_info.mat'));
    
    load(fullfile(sub_paths{i},'intensity_phase'));
    mean_CVC = mean(CVC_phase);
    mean_CVC_smooth = smooth(mean_CVC);
    mean_CVC_smooth_abs = abs(mean_CVC_smooth);
    range1 = [35:75];
    range2 = [76:116];
    range3 = [215:265];
    range4 = [266:360];
    [min_value,pos1] = min(mean_CVC_smooth_abs(range1));
    TA_phase(i) = -180 + (pos1 + 34)+ 15;
    [min_value,pos2] = min(mean_CVC_smooth_abs(range2));
    AP_phase(i) = -180 + (pos2 + 75)+ 15;
    [min_value,pos3] = min(mean_CVC_smooth_abs(range3));
    PD_phase(i) = -180 + (pos3 + 214)+ 15;
    [min_value,pos4] = min(mean_CVC_smooth_abs(range4));
    DT_phase(i) = -180 + (pos4 + 265)+ 15;
    
    
    peak = [AP_phase(i):PD_phase(i)];
    ascend = [TA_phase(i):AP_phase(i)];
    descend = [PD_phase(i):DT_phase(i)];
    trough_bin = ismember(phase_range,trough);
    peak_bin = ismember(phase_range,peak);
    ascend_bin = ismember(phase_range,ascend);
    descend_bin = ismember(phase_range,descend);
    
    GStopo_path = fullfile(sub_paths{i},'fcmat_mean.mat');
    load(GStopo_path)
    fcmat_mean(fcmat_mean < 0) = 0;
    
    trough_mat = mean(fcmat_mean(:,:,trough_bin),3);
    troughmat_R = gretna_inv_fishertrans(trough_mat);
    nanpos = isnan(troughmat_R);
    troughmat_R(nanpos) = 1;
    El_trough(:,i) = efficiency_wei(troughmat_R, 1);
    
    peak_mat = mean(fcmat_mean(:,:,peak_bin),3);
    peakmat_R = gretna_inv_fishertrans(peak_mat);
    nanpos = isnan(peakmat_R);
    peakmat_R(nanpos) = 1;
    El_peak(:,i) = efficiency_wei(peakmat_R, 1);
    
    ascend_mat = mean(fcmat_mean(:,:,ascend_bin),3);
    ascendmat_R = gretna_inv_fishertrans(ascend_mat);
    nanpos = isnan(ascendmat_R);
    ascendmat_R(nanpos) = 1;
    El_ascend(:,i) = efficiency_wei(ascendmat_R, 1);
    
    descend_mat = mean(fcmat_mean(:,:,descend_bin),3);
    descendmat_R = gretna_inv_fishertrans(descend_mat);
    nanpos = isnan(descendmat_R);
    descendmat_R(nanpos) = 1;
    El_descend(:,i) = efficiency_wei(descendmat_R, 1);
end

%% 网络效率 (同样是两个run合起来算的)
clear
clc

sub_paths{1} = 'E:\HCP\Results_task\slow-1';
sub_paths{2} = 'E:\HCP\Results_task\slow-2';
sub_paths{3} = 'E:\HCP\Results_task\slow-3';
sub_paths{4} = 'E:\HCP\Results_task\slow-4';
sub_paths{5} = 'E:\HCP\Results_task\slow-5';

phase_range = [-165:180,-179:-166];

for i = 4:5
    i
    load(fullfile(sub_paths{i},'intensity_phase.mat'));
    mean_CVC = mean(abs(CVC_phase));
    %     mean_CVC = smoothdata(mean_CVC,'gaussian',20);    % 如果直接平均找不到0点，或者有多个0点，就先做一个smooth
    mean_CVC_abs = abs(mean_CVC);
    mean_CVC_abs = [mean_CVC_abs(end-14:end),mean_CVC_abs(1:end-15)];
    range1 = [35:75];
    range2 = [76:116];
    range3 = [215:265];
    range4 = [266:360];
    [min_value,pos1] = min(mean_CVC_abs(range1));
    TA_phase(i) = -180 + (pos1 + 34)+ 15;
    [min_value,pos2] = min(mean_CVC_abs(range2));
    AP_phase(i) = -180 + (pos2 + 75)+ 15;
    [min_value,pos3] = min(mean_CVC_abs(range3));
    PD_phase(i) = -180 + (pos3 + 214)+ 15;
    [min_value,pos4] = min(mean_CVC_abs(range4));
    DT_phase(i) = -180 + (pos4 + 265)+ 15;
    trough = [-179:TA_phase(i),DT_phase(i):180];
    peak = [AP_phase(i):PD_phase(i)];
    ascend = [TA_phase(i):AP_phase(i)];
    descend = [PD_phase(i):DT_phase(i)];
    trough_bin = ismember(phase_range,trough);
    peak_bin = ismember(phase_range,peak);
    ascend_bin = ismember(phase_range,ascend);
    descend_bin = ismember(phase_range,descend);
    
    
    sub_fcmat_path = fullfile(sub_paths{i},'fcmat');
    fcmat_file = dir(sub_fcmat_path);
    fcmat_file(1:2) = [];
    files_name = sort_nat({fcmat_file.name}); % 原函数读取不是按照十进制排序的，这里按照十进制排序，之前写的代码有点问题
    
    for isub = 1:length(fcmat_file)
        load(fullfile(sub_fcmat_path,files_name{isub}));
        fcmat = a_fishertrans(fcmat);   %fisher变换以做不同R值的运算
        trough_mat = mean(fcmat(:,:,trough_bin),3);    %平均波谷
        trough_mat = gretna_inv_fishertrans(trough_mat);     % fisher逆变换以计算网络效率
        nanpos = isnan(trough_mat);
        trough_mat(nanpos) = 0;
        K_value = prctile(reshape(trough_mat,[],1),75);
        trough_mat(trough_mat < K_value) = 0;     %小于0的R值置为0
        El_trough(:,isub) = efficiency_wei(trough_mat, 2);
        Eg_trough(isub) = efficiency_wei(trough_mat);
        
        ascend_mat = mean(fcmat(:,:,ascend_bin),3);    %平均波谷
        ascend_mat = gretna_inv_fishertrans(ascend_mat);     % fisher逆变换以计算网络效率
        nanpos = isnan(ascend_mat);
        ascend_mat(nanpos) = 0;
        K_value = prctile(reshape(ascend_mat,[],1),75);
        ascend_mat(ascend_mat < K_value) = 0;     %小于0的R值置为0
        El_ascend(:,isub) = efficiency_wei(ascend_mat, 2);
        Eg_ascend(isub) = efficiency_wei(ascend_mat);
        
        peak_mat = mean(fcmat(:,:,peak_bin),3);    %平均波谷
        peak_mat = gretna_inv_fishertrans(peak_mat);     % fisher逆变换以计算网络效率
        nanpos = isnan(peak_mat);
        peak_mat(nanpos) = 0;
        K_value = prctile(reshape(peak_mat,[],1),75);
        peak_mat(peak_mat < K_value) = 0;     %小于0的R值置为0
        El_peak(:,isub) = efficiency_wei(peak_mat, 2);
        Eg_peak(isub) = efficiency_wei(peak_mat);
        
        descend_mat = mean(fcmat(:,:,descend_bin),3);    %平均波谷
        descend_mat = gretna_inv_fishertrans(descend_mat);     % fisher逆变换以计算网络效率
        nanpos = isnan(descend_mat);
        descend_mat(nanpos) = 0;
        K_value = prctile(reshape(descend_mat,[],1),75);
        descend_mat(descend_mat < K_value) = 0;     %小于0的R值置为0
        El_descend(:,isub) = efficiency_wei(descend_mat, 2);
        Eg_descend(isub) = efficiency_wei(descend_mat);
    end
    cd(sub_paths{i})
    save Efficiency Eg_trough Eg_ascend Eg_peak Eg_descend El_trough El_ascend El_peak El_descend
end

%% 网络效率和行为表现的相关
clear
clc
sub_paths{1} = 'E:\HCP\Results_task\slow-1';
sub_paths{2} = 'E:\HCP\Results_task\slow-2';
sub_paths{3} = 'E:\HCP\Results_task\slow-3';
sub_paths{4} = 'E:\HCP\Results_task\slow-4';
sub_paths{5} = 'E:\HCP\Results_task\slow-5';


for i = 3:5
    cd(sub_paths{i})
    load('Efficiency.mat');
    load('RT_phase_ACC0.mat')
    load('RT_phase_ACC2.mat')
    T0_ACC = T0;
    T2_ACC = T2;
    A0_ACC = A0;
    A2_ACC = A2;
    P0_ACC = P0;
    P2_ACC = P2;
    D0_ACC = D0;
    D2_ACC = D2;
    load('RT_phase0.mat')
    load('RT_phase2.mat')
    Eg = [Eg_trough,Eg_ascend,Eg_peak,Eg_descend];
    RT0 = [T0(1,:),A0(1,:),P0(1,:),D0(1,:)];
    ACC0 = [T0_ACC(1,:),A0_ACC(1,:),P0_ACC(1,:),D0_ACC(1,:)];
    RT2 = [T2(1,:),A2(1,:),P2(1,:),D2(1,:)];
    ACC2 = [T2_ACC(1,:),A2_ACC(1,:),P2_ACC(1,:),D2_ACC(1,:)];
    [r,p] = corrcoef(Eg,RT0);
    R_global_RT0(i) = r(2);
    P_global_RT0(i) = p(2);
    [r,p] = corrcoef(Eg,ACC0);
    R_global_ACC0(i) = r(2);
    P_global_ACC0(i) = p(2);
    [r,p] = corrcoef(Eg,RT2);
    R_global_RT2(i) = r(2);
    P_global_RT2(i) = p(2);
    [r,p] = corrcoef(Eg,ACC2);
    R_global_ACC2(i) = r(2);
    P_global_ACC2(i) = p(2);
    
    cd(sub_paths{i})
    save RT_ACC_Eg RT0 RT2 ACC0 ACC2 Eg

    
    for iROI = 1:246
        El = [El_trough(iROI,:),El_ascend(iROI,:),El_peak(iROI,:),El_descend(iROI,:)];
        [r,p] = corrcoef(El,RT0);
        R_local_RT0(iROI,i) = r(2);
        p_local_RT0(iROI,i) = p(2);
        [r,p] = corrcoef(El,ACC0);
        R_local_ACC0(iROI,i) = r(2);
        p_local_ACC0(iROI,i) = p(2);
        
        [r,p] = corrcoef(El,RT2);
        R_local_RT2(iROI,i) = r(2);
        p_local_RT2(iROI,i) = p(2);
        [r,p] = corrcoef(El,ACC2);
        R_local_ACC2(iROI,i) = r(2);
        p_local_ACC2(iROI,i) = p(2);
    end
    
    mode = 1;
    FDR_local_ACC0(:,i) = a_multicorrect(R_local_ACC0(:,i),p_local_ACC0(:,i),mode);
    FDR_local_RT0(:,i) = a_multicorrect(R_local_RT0(:,i),p_local_RT0(:,i),mode);
    FDR_local_ACC2(:,i) = a_multicorrect(R_local_ACC2(:,i),p_local_ACC2(:,i),mode);
    FDR_local_RT2(:,i) = a_multicorrect(R_local_RT2(:,i),p_local_RT2(:,i),mode);
end
cd('E:\HCP\Results_task')
save FDR_El_RT_ACC FDR_local_RT0 FDR_local_RT2 FDR_local_ACC0 FDR_local_ACC2

A = mean(Eg_trough);
B = mean(Eg_ascend);
C = mean(Eg_peak);
D = mean(Eg_descend);

Eg = [Eg_trough,Eg_ascend,Eg_peak,Eg_descend];
RT = [T0_ACC(1,:),A0_ACC(1,:),P0_ACC(1,:),D0_ACC(1,:)];
[r_t,p_t] = corrcoef(Eg,RT);

%% task vs. rest
clear
clc

sub_paths{1} = 'E:\HCP\Results_PL\slow-1';
sub_paths{2} = 'E:\HCP\Results_PL\slow-2';
sub_paths{3} = 'E:\HCP\Results_PL\slow-3';
sub_paths{4} = 'E:\HCP\Results_PL\slow-4';
sub_paths{5} = 'E:\HCP\Results_PL\slow-5';

sub_paths2{1} = 'E:\HCP\Results_task\slow-1';
sub_paths2{2} = 'E:\HCP\Results_task\slow-2';
sub_paths2{3} = 'E:\HCP\Results_task\slow-3';
sub_paths2{4} = 'E:\HCP\Results_task\slow-4';
sub_paths2{5} = 'E:\HCP\Results_task\slow-5';

for i = 1:5
    load(fullfile(sub_paths{i},'intensity_phase.mat'));
    rest_CVC = CVC_phase;
    for icl = 1:360
        x = rest_CVC(:,icl);
        x(icl) = [];
        mean_cl_z1(icl,:) = x;
    end
    load(fullfile(sub_paths2{i},'intensity_phase.mat'));
    task_CVC = CVC_phase;
    for icl = 1:360
        x = task_CVC(:,icl);
        x(icl) = [];
        mean_cl_z2(icl,:) = x;
    end
    [r,p] = corrcoef(mean_cl_z1,mean_cl_z2);
    R(i) = r(2);
    P(i) = p(2);
end
