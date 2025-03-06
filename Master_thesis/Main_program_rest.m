% 重采样呼吸、心跳数据到fMRI扫描频率，再线性回归每个体素时间序列的呼吸、心跳、（代码）
% 6个头动参数和它们的一阶导数、以及每次扫描的线性趋势（代码）
% 白质、脑脊液信号（dpabi）

%%
% 用dpabi做平滑

%% FIR滤波
clear

bands = ccs_core_lfobands(1200, 0.72); % 分频段
% for iband = 1:4
%     bands{iband} = roundn(bands{iband},-4);
% end
% bands{1} = [];
% bands{2} = [];
% bands{3} = [];
% bands(cellfun(@isempty,bands))=[];
% bands{1} = [0.2225,0.24];

sub_path = 'F:\lifespan\FunRawARWSC';
out_path = 'F:\lifespan\FunRawARWSCF';
Fs = 1/2;

a_Nii_FIRfilter(sub_path,out_path,bands,Fs); % 滤波函数

%% 提取ROI，GS信号
clear
clc
sub_paths{1} = 'E:\HCP\Results_PL\slow-1';
sub_paths{2} = 'E:\HCP\Results_PL\slow-2';
sub_paths{3} = 'E:\HCP\Results_PL\slow-3';
sub_paths{4} = 'E:\HCP\Results_PL\slow-4';
sub_paths{5} = 'E:\HCP\Results_PL\slow-5';

out_paths{1} = 'E:\HCP\Results_PL\slow-1';
out_paths{2} = 'E:\HCP\Results_PL\slow-2';
out_paths{3} = 'E:\HCP\Results_PL\slow-3';
out_paths{4} = 'E:\HCP\Results_PL\slow-4';
out_paths{5} = 'E:\HCP\Results_PL\slow-5';

mask_path = 'D:\program\matlab\Common\Mask\BN_Atlas_246_2mm.nii';

for i = 1:5
    [ROISignals,GS] = a_GSROI(sub_paths{i},mask_path);
    mkdir(out_paths{i})
    cd(out_paths{i});
    save ROISignals ROISignals
    save GS GS
end

%% 计算每个相位上的GS拓扑图
clear
clc
sub_paths{1} = 'E:\HCP\Results_PL\slow-1';
sub_paths{2} = 'E:\HCP\Results_PL\slow-2';
sub_paths{3} = 'E:\HCP\Results_PL\slow-3';
sub_paths{4} = 'E:\HCP\Results_PL\slow-4';
sub_paths{5} = 'E:\HCP\Results_PL\slow-5';

out_paths{1} = 'E:\HCP\Results_PL\slow-1';
out_paths{2} = 'E:\HCP\Results_PL\slow-2';
out_paths{3} = 'E:\HCP\Results_PL\slow-3';
out_paths{4} = 'E:\HCP\Results_PL\slow-4';
out_paths{5} = 'E:\HCP\Results_PL\slow-5';
% GStopo and FCmat
for i = 1:5
    i
    load(fullfile(sub_paths{i},'GS.mat'))
    % 先计算GS相位
    GS_phase = a_Phase(GS);
    load(fullfile(sub_paths{i},'GS_phase.mat'))
    load(fullfile(sub_paths{i},'ROISignals.mat'))

    sldwd_GStopo(GS,GS_phase,ROISignals,5,0.5,out_paths{i},1) % 计算每个相位上的GS拓扑图
    sldwd_fcmat(GS,GS_phase,ROISignals,5,0.5,out_paths{i},1) % 计算每个相位上的两两脑区的FC
end

%% 计算平均GS拓扑图和FCmat 
clear
clc
sub_paths{1} = 'E:\HCP\Results_PL\slow-1';
sub_paths{2} = 'E:\HCP\Results_PL\slow-2';
sub_paths{3} = 'E:\HCP\Results_PL\slow-3';
sub_paths{4} = 'E:\HCP\Results_PL\slow-4';
sub_paths{5} = 'E:\HCP\Results_PL\slow-5';

out_paths{1} = 'E:\HCP\Results_PL\slow-1';
out_paths{2} = 'E:\HCP\Results_PL\slow-2';
out_paths{3} = 'E:\HCP\Results_PL\slow-3';
out_paths{4} = 'E:\HCP\Results_PL\slow-4';
out_paths{5} = 'E:\HCP\Results_PL\slow-5';


% GStopo intensity and correlation
for i = 1:5
    i
    cd(fullfile(sub_paths{i},'sub_topo'));
    sub_file = dir(fullfile(sub_paths{i},'sub_topo'));
    sub_file(1:2) = [];
    GS_topoZ_allsubj = zeros(246,360,328);
    for isubj = 1:328
        load(sub_file(isubj).name)
        GS_topoZ = a_fishertrans(GS_topo);
        GS_topoZ_allsubj(:,:,isubj) = GS_topoZ;
    end
    GS_topo_mean = mean(GS_topoZ_allsubj,3);
    cd(out_paths{i})
    save GS_topo_mean GS_topo_mean
end

% fcmat intensity and correlation
for i = 1:5
    i
    cd(fullfile(sub_paths{i},'fcmat'));
    sub_file = dir(fullfile(sub_paths{i},'fcmat'));
    sub_file(1:2) = [];
    fcmatZ_allsubj = zeros(246,246,360,328);
    for isubj = 1:328
        isubj
        load(sub_file(isubj).name)
        fcmatZ = a_fishertrans(fcmat);
        fcmatZ_allsubj(:,:,:,isubj) = fcmatZ;
    end
    fcmat_mean = mean(fcmatZ_allsubj,4);
    cd(out_paths{i})
    save fcmat_mean fcmat_mean
end

%% display strength and correlation
clear
clc
sub_paths{1} = 'E:\HCP\Results_PL\slow-1';
sub_paths{2} = 'E:\HCP\Results_PL\slow-2';
sub_paths{3} = 'E:\HCP\Results_PL\slow-3';
sub_paths{4} = 'E:\HCP\Results_PL\slow-4';
sub_paths{5} = 'E:\HCP\Results_PL\slow-5';

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

sub_paths{1} = 'E:\HCP\Results_PL\slow-1';
sub_paths{2} = 'E:\HCP\Results_PL\slow-2';
sub_paths{3} = 'E:\HCP\Results_PL\slow-3';
sub_paths{4} = 'E:\HCP\Results_PL\slow-4';
sub_paths{5} = 'E:\HCP\Results_PL\slow-5';

for i = 1:5
    load(fullfile(sub_paths{i},'GS_topo_mean.mat'))
    GS_topo_mean = [GS_topo_mean(:,end-14:end),GS_topo_mean(:,1:end-15)];
    GS_topo_mean = GS_topo_mean';
    GS_topo_mean_hilb = hilbert(GS_topo_mean);
    strength_phase = angle(GS_topo_mean_hilb);
    strength_phase = strength_phase * 180 / pi; % 算相位
    
    GS_topo_mean_all = mean(GS_topo_mean,2);
    GS_topo_mean_all_hilb = hilbert(GS_topo_mean_all);
    strength_phase_all = angle(GS_topo_mean_all_hilb);
    strength_phase_all = strength_phase_all * 180 / pi;
    
    phase_difference = bsxfun(@minus,strength_phase,strength_phase_all); % 减去平均相位
    CVC_phase = corr(phase_difference',phase_difference','type','pearson'); % 求两两相位的相似性
    cd(sub_paths{i});
    save intensity_phase phase_difference CVC_phase
end

%% 区分4种状态：波峰，波谷，上升，下降
clear
clc

sub_paths{1} = 'E:\HCP\Results_PL\slow-1';
sub_paths{2} = 'E:\HCP\Results_PL\slow-2';
sub_paths{3} = 'E:\HCP\Results_PL\slow-3';
sub_paths{4} = 'E:\HCP\Results_PL\slow-4';
sub_paths{5} = 'E:\HCP\Results_PL\slow-5';


for i = 1:5
    i
    load(fullfile(sub_paths{i},'intensity_phase.mat'));
    for icl = 1:360
        x = CVC_phase(:,icl);
        x(icl) = [];
        mean_cl_z(icl) = mean(abs(a_fishertrans(x)));
    end
    mean_cl_r = gretna_inv_fishertrans(mean_cl_z);
    cd(sub_paths{i})
    
    mean_cl_p = a_r2p(mean_cl_r,246).*360; % 使r值变P值 FWE校正
    x = mean_cl_r;
    x(mean_cl_p > 0.05) = 1;
    min_r(i) = min(abs(x));
    %     mean_CVC = smoothdata(mean_CVC,'gaussian',20);    % 如果直接平均找不到0点，或者有多个0点，就先做一个smooth
    save mean_cl_r mean_cl_r
end
    
    


%% 网络效率
clear
clc

sub_paths{1} = 'E:\HCP\Results_PL\slow-1';
sub_paths{2} = 'E:\HCP\Results_PL\slow-2';
sub_paths{3} = 'E:\HCP\Results_PL\slow-3';
sub_paths{4} = 'E:\HCP\Results_PL\slow-4';
sub_paths{5} = 'E:\HCP\Results_PL\slow-5';

phase_range = [-180:179];

for i = 1:5
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
    
    
    range1 = [35:75];
    range2 = [76:120];
    range3 = [215:265];
    range4 = [266:360];
    [min_value,pos1] = min(mean_cl_r(range1));
    TA_phase(i) = -180 + (pos1 + 34);
    [min_value,pos2] = min(mean_cl_r(range2));
    AP_phase(i) = -180 + (pos2 + 75);
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

% 计算全局效率和局部效率
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

sub_paths{1} = 'E:\HCP\Results_PL\slow-1';
sub_paths{2} = 'E:\HCP\Results_PL\slow-2';
sub_paths{3} = 'E:\HCP\Results_PL\slow-3';
sub_paths{4} = 'E:\HCP\Results_PL\slow-4';
sub_paths{5} = 'E:\HCP\Results_PL\slow-5';

% 被试内方差分析
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

eta_sq = a_etasq(F,3,981); % 效应量

% 局部效率的方差分析
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
    
%     mean_El_ascend(:,i) = mean(El_ascend,2);
%     mean_El_peak(:,i) = mean(El_peak,2);
%     mean_El_descend(:,i) = mean(El_descend,2);
%     mean_El_trough(:,i) = mean(El_trough,2);
%     
%     for iROI = 1:246
%         [h,p,ci,stats] = ttest(El_ascend(iROI,:));
%         El_ascend_T(iROI,i) = stats.tstat;
%         [h,p,ci,stats] = ttest(El_peak(iROI,:));
%         El_peak_T(iROI,i) = stats.tstat;
%         [h,p,ci,stats] = ttest(El_descend(iROI,:));
%         El_descend_T(iROI,i) = stats.tstat;
%         [h,p,ci,stats] = ttest(El_trough(iROI,:));
%         El_trough_T(iROI,i) = stats.tstat;
%     end
    
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








