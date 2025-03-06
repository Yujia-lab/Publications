%% cluster dimensionality reduction
clear
clc
cohere_path = 'E:\GS Coherence\result\data\GS_topo_Cohere.mat';
load(cohere_path);
age_path = 'E:\GS Coherence\result\data\AGECORR.mat';
load(age_path);


mean_cohere = mean(ROI_Cohere,3);
label1 = find(FDR>0);
label2 = find(FDR<0);

FDR2=FDR;
FDR2(label1) = 1;
FDR2(label2) = -1;

opts = statset('MaxIter', 500, 'Display', 'off');
[IDX, C, SUMD, D] = kmeans(FDR2', 3,'Replicates',5,'options',opts);

obj = evalclusters(FDR2', 'kmeans', 'DaviesBouldin','klist',2:8);
obj = evalclusters(Z', 'kmeans', 'DaviesBouldin','klist',2:20);

obj = evalclusters(mean_cohere, 'kmeans', 'CalinskiHarabasz','klist',1:6);
obj = evalclusters(Z, 'kmeans', 'CalinskiHarabasz','klist',1:8);
[IDX, C, SUMD, D] = kmeans(mean_cohere, 4);



mappedX = tsne(mean_cohere);
C1 = find(IDX == 1);
C2 = find(IDX == 2);
C3 = find(IDX == 3);
C4 = find(IDX == 4);
scatter(mappedX(C1,1), mappedX(C1,2));
hold on 
scatter(mappedX(C2,1), mappedX(C2,2));
scatter(mappedX(C3,1), mappedX(C3,2));
scatter(mappedX(C4,1), mappedX(C4,2));

cohere_path = 'D:\Aoyujia\Academic\Writing\Lifespan Coherence\GS_topo_Cohere.mat';
load(cohere_path);

mean_cohere = mean(ROI_Cohere,3);
[IDX, C, SUMD, D] = kmeans(mean_cohere', 2);


mappedX = tsne(mean_cohere');
C1 = find(IDX == 1);
C2 = find(IDX == 2);
% C3 = find(IDX == 3);
% C4 = find(IDX == 4);
scatter(mappedX(C1,1), mappedX(C1,2));
hold on 
scatter(mappedX(C2,1), mappedX(C2,2));
% scatter(mappedX(C3,1), mappedX(C3,2));
% scatter(mappedX(C4,1), mappedX(C4,2));

[r,p] = corrcoef(C(1,50:100),C(2,50:100));

%% vector 2 nii
clc
clear
cohere_path = 'E:\GS Coherence\result\allsubcohere\GS_topo_Cohere.mat';
load(cohere_path);
cluster_path = 'E:\GS Coherence\result\cluster2.0\CORR_freq\cluster_frequency.mat';
load(cluster_path);
C1 = 14;
C2 = 63;
C3 = 160;
mask_path = 'D:\program\matlab\Common\Mask\BN_Atlas_246_3mm.nii';
mask_data = y_Read(mask_path);
brain_path = 'F:\lifespan\FunRawARWSC\sub-031275\CovRegressed_4DVolume.nii';
[brain,header] = y_Read(brain_path,1);



for i = 1:322
    if i < 10
        name = ['00',num2str(i)];
    elseif i >= 10 & i < 100
        name = ['0',num2str(i)];   
    else
        name = num2str(i);
    end
    brain1 = mean(ROI_Cohere(1:C1,:,i),1);
    out_path1 = 'E:\lifespan_phaselag\fre1';
    mkdir(out_path1)
    cd(out_path1)
    a_vec2nii(brain1,mask_data,header,name)
    
    brain2 = mean(ROI_Cohere(1+C1:C2,:,i),1);
    out_path2 = 'E:\lifespan_phaselag\fre2';
    mkdir(out_path2)
    cd(out_path2)
    a_vec2nii(brain2,mask_data,header,name)
    
    brain3 = mean(ROI_Cohere(1+C2:C3,:,i),1);
    out_path3 = 'E:\lifespan_phaselag\fre3';
    mkdir(out_path3)
    cd(out_path3)
    a_vec2nii(brain3,mask_data,header,name)
    
    brain4 = mean(ROI_Cohere(1+C3:end,:,i),1);
    out_path4 = 'E:\lifespan_phaselag\fre4';
    mkdir(out_path4)
    cd(out_path4)
    a_vec2nii(brain4,mask_data,header,name)
end
    
%% 平均Coherence算的相关
clc
clear
cohere_path = 'E:\GS Coherence\result\allsubcohere\GS_topo_Cohere.mat';
load(cohere_path);
cluster_path = 'E:\GS Coherence\result\cluster2.0\CORR_freq\cluster_frequency.mat';
load(cluster_path);
info_path = 'F:\lifespan\sub_information_deleted.xlsx';
info = xlsread(info_path);
age=info(:,3);
gender=info(:,4); 

C1 = 14;
C2 = 63;
C3 = 160;
mask_path = 'D:\program\matlab\Common\Mask\BN_Atlas_246_3mm.nii';
mask_data = y_Read(mask_path);
brain_path = 'F:\lifespan\FunRawARWSC\sub-031275\CovRegressed_4DVolume.nii';
[brain,header] = y_Read(brain_path,1);

X=[gender,ones(length(info),1)];  

for iROI = 1:246
    ROI_Cohere2 = squeeze(ROI_Cohere(:,iROI,:));
    mean_Cohere1 = mean(ROI_Cohere2(1:C1,:),1);
    mean_Cohere2 = mean(ROI_Cohere2(C1+1:C2,:),1);
    mean_Cohere3 = mean(ROI_Cohere2(C2+1:C3,:),1);
    mean_Cohere4 = mean(ROI_Cohere2(C3+1:end,:),1);
    [b,bint,Re1,Rint,stats]=regress(mean_Cohere1',X);
    [b,bint,Re2,Rint,stats]=regress(mean_Cohere2',X);
    [b,bint,Re3,Rint,stats]=regress(mean_Cohere3',X);
    [b,bint,Re4,Rint,stats]=regress(mean_Cohere4',X);
    [r1,p1]=corrcoef(Re1,age);
    [r2,p2]=corrcoef(Re2,age);
    [r3,p3]=corrcoef(Re3,age);
    [r4,p4]=corrcoef(Re4,age);
    R(iROI,1) = r1(2);
    R(iROI,2) = r2(2);
    R(iROI,3) = r3(2);
    R(iROI,4) = r4(2);
    P(iROI,1) = p1(2);
    P(iROI,2) = p2(2);
    P(iROI,3) = p3(2);
    P(iROI,4) = p4(2);
end

out_path = 'E:\GS Coherence\result\age-related_GStopo';
for i = 1:4
    FDR(i,:) = a_FDR(R(:,i),P(:,i));
    cd(out_path);
    a_vec2nii(FDR(i,:),mask_data,header,['FR',num2str(i)]);
end

%% for neurosynth
clear
clc
path = 'E:\GS Coherence\result\cluster2.0\CORR_freq\cluster_frequency.mat';
load(path);
out_path = 'E:\GS Coherence\result\cluster2.0\CORR_freq';
cd(out_path)
for i = 1:4
    a_vec2nii(mean_data(i,:),mask_data,header,['unthreshold',num2str(i)]);
end
%% correlation between GStopo and agetopo
clc
clear
cohere_path = 'E:\GS Coherence\result\GS_topo_Cohere.mat';
load(cohere_path);
age_path = 'E:\GS Coherence\result\data\AGECORR.mat';
load(age_path);
mask_path = 'D:\program\matlab\Common\Mask\BN_Atlas_246_3mm.nii';
mask_data = y_Read(mask_path);
brain_path = 'F:\lifespan\FunRawARWSC\sub-031275\CovRegressed_4DVolume.nii';
[brain,header] = y_Read(brain_path,1);

cd('E:\GS Coherence\result');

x = reshape(ROI_Cohere,257*246,1);
y = reshape(Z',257*246,1);
[r,p] = corrcoef(x,y);
R_all = r(2);
P_all = p(2);
%时间去分化
for i = 1:246
    [r,p] = corrcoef(ROI_Cohere(i,:),Z(:,i));
    R_temp(i) = r(2);
    P_temp(i) = p(2);
end

%空间去分化
for i = 1:257
    [r,p] = corrcoef(ROI_Cohere(:,i),Z(i,:));
    R_betw(i) = r(2);
    P_betw(i) = p(2);
end

% 分别计算时间去分化的正负R值
R_tempfdr = a_FDR(R_temp,P_temp);
Z = Z';
C_pos = ROI_Cohere(R_tempfdr>0,:);
C_neg = ROI_Cohere(R_tempfdr<0,:);
x2 = reshape(C_pos,257*59,[]);
x3 = reshape(C_neg,257*164,[]);
R_pos = Z(R_tempfdr>0,:);
R_neg = Z(R_tempfdr<0,:);
y2 = reshape(R_pos,257*59,[]);
y3 = reshape(R_neg,257*164,[]);

[r,p] = corrcoef(x2,y2);
R_pos = r(2);
P_pos = p(2);
[r,p] = corrcoef(x3,y3);
R_neg = r(2);
P_neg = p(2);


% 分别计算空间去分化的正负R值
FDR_r = a_FDR(R_betw,P_betw);
C1 = mean(ROI_Cohere(:,1:8),2);
Z1 = mean(Z(:,1:8),2);
C2 = mean(ROI_Cohere(:,27:46),2);
Z2 = mean(Z(:,27:46),2);
C3 = mean(ROI_Cohere(:,164:end),2);
Z3 = mean(Z(:,164:end),2);
[r,p] = corrcoef(C1,Z1);
R1 = r(2);
P1 = p(2);
[r,p] = corrcoef(C2,Z2);
R2 = r(2);
P2 = p(2);
[r,p] = corrcoef(C3,Z3);
R3 = r(2);
P3 = p(2);


cd('E:\GS Coherence\result\cluster2.0\Dediff')
save spatial_dediff.mat R_temp P_temp R_betw P_betw R_tempfdr FDR_r
%% 4 clusters CVC
clc
clear

CORR_path = 'E:\lifespan_phaselag\result\cluster2.0\CORR_freq\cluster_frequency.mat';
Cohere_path = 'E:\lifespan_phaselag\result\cluster2.0\Cohere_freq\cluster_frequency.mat';

load(Cohere_path);
mean_data1 = mean_data;
load(CORR_path);
mean_data2 = mean_data;
for i1 = 1:4
    for i2 = 1:4
        [r,p] = corrcoef(mean_data1(i1,:),mean_data1(i2,:));
        R1(i1,i2) = r(2);
        [r,p] = corrcoef(mean_data2(i1,:),mean_data2(i2,:));
        R2(i1,i2) = r(2);
    end
end

%% 结果报告
clear
clc
load('E:\GS Coherence\result\cluster2.0\CORR_freq\cluster_frequency.mat');

for i = 1:246
    FR1 = FDR_Z(1,i);
    FR3 = FDR_Z(3,i);
    if FR1 == FR3
        diff_mark(i) = 0;
    elseif FR1 < 0 & FR3 < 0 
            diff_mark(i) = -2;
        elseif FR1 < 0 & FR3 == 0
            diff_mark(i) = -1;
            elseif FR1 == 0 & FR3 < 0
            diff_mark(i) = 1;
    else
        diff_mark(i) = 0;
    end
end
        
%% 计算肯德尔和谐系数
clc
clear
load('E:\GS Coherence\result\data\GS_topo_Cohere.mat');
load('E:\GS Coherence\result\cluster2.0\Cohere_freq\cluster_frequency.mat');
FR1 = 23;
FR2 = 68;
FR3 = 121;
FR4 = 257;

mean_Cohere = mean(ROI_Cohere,3)';
Rank_Cohere = order2grade(mean_Cohere);

Rank_mean = mean(Rank_Cohere, 2);
Rank_base = bsxfun(@minus,Rank_Cohere,Rank_mean);
obj = evalclusters(Rank_Cohere', 'kmeans', 'CalinskiHarabasz','klist',1:8);
opts = statset('MaxIter', 500, 'Display', 'off');
[IDX, C, SUMD, D] = kmeans(Rank_Cohere', 4,'Replicates',5,'options',opts);

cd('E:\GS Coherence\result\cluster2.0\Cohere_freq')
save rank_cluster.mat IDX C 

a_vec2nii(IDX,'Rank_cluster.nii');



% FR1_data = squeeze(mean(ROI_Cohere(1:FR1,:,:),1));
% FR2_data = squeeze(mean(ROI_Cohere(FR1+1:FR2,:,:),1));
% FR3_data = squeeze(mean(ROI_Cohere(FR2+1:FR3,:,:),1));
% FR4_data = squeeze(mean(ROI_Cohere(FR3+1:FR4,:,:),1));
% 
% Rank_FR1 = order2grade(FR1_data);
% Rank_FR2 = order2grade(FR2_data);
% Rank_FR3 = order2grade(FR3_data);
% Rank_FR4 = order2grade(FR4_data);
% 
% for i = 1:246
%     ROI_FR1 = Rank_FR1(i,:);
%     ROI_FR2 = Rank_FR2(i,:);
%     ROI_FR3 = Rank_FR3(i,:);
%     ROI_FR4 = Rank_FR4(i,:);
%     all_FR = [ROI_FR1',ROI_FR2',ROI_FR3',ROI_FR4'];
%     [p,tbl,stat] = friedman(all_FR,1,'off');
%     P_friedman(i) = p;
%     chi(i) = tbl{2,5};     
%     w(i) = kandallw(all_FR);
% end
%
% FDR = mafdr(P_friedman, 'BHFDR',true);
% neg = find(FDR > 0.05);
% select_ROI = find(FDR < 0.05);
% 
% FR1_mean = squeeze(mean(mean_Cohere(1:FR1,:),1));
% FR2_mean = squeeze(mean(mean_Cohere(FR1+1:FR2,:),1));
% FR3_mean = squeeze(mean(mean_Cohere(FR2+1:FR3,:),1));
% FR4_mean = squeeze(mean(mean_Cohere(FR3+1:FR4,:),1));
% all_mean = [FR1_mean',FR2_mean',FR3_mean',FR4_mean'];
% Rank_mean = order2grade(all_mean);
% 
% for i = 1:3
%     Rank_change(:,i) = Rank_mean(:,i+1) - Rank_mean(:,i);
% end
%      
% obj = evalclusters(Rank_change(select_ROI,:), 'kmeans', 'CalinskiHarabasz','klist',1:10);
% 
% opts = statset('MaxIter', 500, 'Display', 'off');
% [IDX, C, SUMD, D] = kmeans(Rank_change(select_ROI,:), 2,'Replicates',5,'options',opts);


%% 画皮层下的图
clear
clc
cmap_path = 'E:\GS Coherence\Colormap\coolwarm.csv';
ColorMap = csvread(cmap_path);
cd('E:\GS Coherence\figure-sub')
save coolwarm ColorMap

Cohere_path = 'E:\GS Coherence\result\GS_topo_Cohere.mat';
load(Cohere_path);
mean_cohere = mean(ROI_Cohere,2);
mean_cohere(1:210) = 0;
cd('E:\GS Coherence\result\cluster2.0\Cohere_freq')
a_vec2nii(mean_cohere,'Cohereall_sub');

Cohere_path = 'E:\GS Coherence\result\cluster2.0\Cohere_freq\cluster_frequency.mat';
load(Cohere_path);
mean_data(:,1:210) = 0;
cd('E:\GS Coherence\result\cluster2.0\Cohere_freq')
a_vec2nii(mean_data(1,:),'Cohere4-1_sub');
a_vec2nii(mean_data(2,:),'Cohere4-2_sub');
a_vec2nii(mean_data(3,:),'Cohere4-3_sub');
a_vec2nii(mean_data(4,:),'Cohere4-4_sub');

% rank_path = 'E:\GS Coherence\result\cluster2.0\Cohere_freq\rank_cluster.mat';
% load(rank_path);
% rank_cluster = IDX;
% rank_cluster(IDX == 1) = -0.9;
% rank_cluster(IDX == 2) = 0.9;
% rank_cluster(1:210) = 0;
% cd('E:\GS Coherence\result\cluster2.0\Cohere_freq')
% a_vec2nii(rank_cluster,'rank_sub');


load('E:\GS Coherence\result\cluster2.0\CORR_freq\CORR_all.mat')
FDR_Z(:,1:210) = 0;
cd('E:\GS Coherence\result\cluster2.0\CORR_freq')
a_vec2nii(FDR_Z,'CORRall_sub');

CORR_path = 'E:\GS Coherence\result\cluster2.0\CORR_freq\cluster_frequency.mat';
load(CORR_path);
FDR_Z(:,1:210) = 0;
cd('E:\GS Coherence\result\cluster2.0\CORR_freq')
a_vec2nii(FDR_Z(1,:),'CORR4-1_sub');
a_vec2nii(FDR_Z(2,:),'CORR4-2_sub');
a_vec2nii(FDR_Z(3,:),'CORR4-3_sub');
a_vec2nii(FDR_Z(4,:),'CORR4-4_sub');

betw_path = 'E:\GS Coherence\result\data\Corr_betw.mat';
load(betw_path);
CORR_betw_FDR(1:210) = 0;
cd('E:\GS Coherence\result')
a_vec2nii(CORR_betw_FDR,'betw_sub');

%% 画两个cluster的图
path = 'E:\GS Coherence\result\cluster2.0\Cohere_freq\rank_cluster.mat';
load(path);
C = C';
g=gramm('x',X3,'y',Y3);

%% 全频段的CORR
load('E:\GS Coherence\result\data\AGECORR.mat')
mean_Z = mean(Z,1);
mean_R = gretna_inv_fishertrans(mean_Z);
mean_P = a_r2p(mean_R,322);
FDR_Z = a_FDR(mean_Z,mean_P);

cd('E:\GS Coherence\result\cluster2.0\CORR_freq')
a_vec2nii(mean_Z,'unthreshold_all')


cd('E:\GS Coherence\result\cluster2.0\CORR_freq')
save CORR_all.mat FDR_Z

%% 强度归一化做ANOVA
clear
path = 'E:\GS Coherence\result\data\GS_topo_Cohere.mat';
load(path)

Cohere_fre(1,:,:) = mean(ROI_Cohere(1:23,:,:),1);
Cohere_fre(2,:,:) = mean(ROI_Cohere(24:45,:,:),1);
Cohere_fre(3,:,:) = mean(ROI_Cohere(46:136,:,:),1);
Cohere_fre(4,:,:) = mean(ROI_Cohere(137:end,:,:),1);
Cohere_fre(5,:,:) = mean(ROI_Cohere,1);

for isub = 1:322
    for ifre = 1:5
        Cohere_fre1 = Cohere_fre(ifre,:,isub);
        av = mean(Cohere_fre1);
        Cohere_norm(ifre,:,isub) = Cohere_fre1 - av;
    end
end

group(1:322,1) = zeros(322,1);
group = [group;ones(322,1)];
group = [group;2*ones(322,1)];
group = [group;3*ones(322,1)];
for iROI = 1:246
    Cohere_ROI = squeeze(Cohere_norm(1,iROI,:));
    Cohere_ROI = [Cohere_ROI;squeeze(Cohere_norm(2,iROI,:))];
    Cohere_ROI = [Cohere_ROI;squeeze(Cohere_norm(3,iROI,:))];
    Cohere_ROI = [Cohere_ROI;squeeze(Cohere_norm(4,iROI,:))];
    [d, p, stats] = manova1(Cohere_ROI,group);
    P(iROI) = p;
    F(iROI) = (stats.B/stats.dfB)/(stats.W/stats.dfW);
end
F_FDR = a_FDR(F,P);

ETAsq = F_FDR.*stats.dfB./(F_FDR.*stats.dfB+stats.dfW);
cohenf = sqrt(ETAsq./(1-ETAsq));



for iROI = 1:246
    [h,p1,ci,stats1] = ttest(Cohere_fre(1,iROI,:),Cohere_fre(2,iROI,:));
    [h,p2,ci,stats2] = ttest(Cohere_fre(1,iROI,:),Cohere_fre(3,iROI,:));
    [h,p3,ci,stats3] = ttest(Cohere_fre(1,iROI,:),Cohere_fre(4,iROI,:));
    [h,p4,ci,stats4] = ttest(Cohere_fre(2,iROI,:),Cohere_fre(3,iROI,:));
    [h,p5,ci,stats5] = ttest(Cohere_fre(2,iROI,:),Cohere_fre(4,iROI,:));
    [h,p6,ci,stats6] = ttest(Cohere_fre(3,iROI,:),Cohere_fre(4,iROI,:));
    P(1,iROI) = p1;
    P(2,iROI) = p2;
    P(3,iROI) = p3;
    P(4,iROI) = p4;
    P(5,iROI) = p5;
    P(6,iROI) = p6;
    T(1,iROI) = stats1.tstat;
    T(2,iROI) = stats2.tstat;
    T(3,iROI) = stats3.tstat;
    T(4,iROI) = stats4.tstat;
    T(5,iROI) = stats5.tstat;
    T(6,iROI) = stats6.tstat;
    ES1(1,iROI) = stats1.tstat/sqrt(321);
    ES2(1,iROI) = stats2.tstat/sqrt(321);
    ES3(1,iROI) = stats3.tstat/sqrt(321);
    ES4(1,iROI) = stats4.tstat/sqrt(321);
    ES5(1,iROI) = stats5.tstat/sqrt(321);
    ES6(1,iROI) = stats6.tstat/sqrt(321);
end
ES = [ES1;ES2;ES3;ES4;ES5;ES6];
FDR(1,:) = mafdr(P(1,:),'BHFDR',true);
FDR(2,:) = mafdr(P(2,:),'BHFDR',true);
FDR(3,:) = mafdr(P(3,:),'BHFDR',true);
FDR(4,:) = mafdr(P(4,:),'BHFDR',true);
FDR(5,:) = mafdr(P(5,:),'BHFDR',true);
FDR(6,:) = mafdr(P(6,:),'BHFDR',true);
sig_pos = find(FDR>0.05);
ES(sig_pos) = 0;
cd('E:\GS Coherence\result\cluster2.0\Cohere_freq')
save ES.mat ES

% 画皮下的图
ES(:,1:210) = 0;
cd('E:\GS Coherence\result\cluster2.0\Cohere_freq')
a_vec2nii(ES(1,:),'ES1');
a_vec2nii(ES(2,:),'ES2');
a_vec2nii(ES(3,:),'ES3');
a_vec2nii(ES(4,:),'ES4');
a_vec2nii(ES(5,:),'ES5');
a_vec2nii(ES(6,:),'ES6');




for i = 1:6
    for ii = 1:6
    [r,p] = corrcoef(T(i,:),T(ii,:));
    R(i,ii) = r(2);
    end
end

% 校正
FDR(1,:) = mafdr(P(1,:),'BHFDR',true);
FDR(2,:) = mafdr(P(2,:),'BHFDR',true);
FDR(3,:) = mafdr(P(3,:),'BHFDR',true);
FDR(4,:) = mafdr(P(4,:),'BHFDR',true);
FDR(5,:) = mafdr(P(5,:),'BHFDR',true);
FDR(6,:) = mafdr(P(6,:),'BHFDR',true);

sig_pos = find(FDR>0.05);
T(sig_pos) = 0;
%二值化方便查看
%T(T<0) = -1;
%T(T>0) = 1;
T(:,[1:210]) = 0;





brain0 = zeros(1,246);
brain0(rank4) = 1;
cd('E:\GS Coherence\result\cluster2.0\Cohere_freq')
a_vec2nii(brain0,'rank4');

brain0 = zeros(1,246);
brain0(rank3) = 1;
cd('E:\GS Coherence\result\cluster2.0\Cohere_freq')
a_vec2nii(brain0,'rank3');

brain0 = zeros(1,246);
brain0(rank2) = 1;
cd('E:\GS Coherence\result\cluster2.0\Cohere_freq')
a_vec2nii(brain0,'rank2');

brain0 = zeros(1,246);
brain0(rank1) = 1;
cd('E:\GS Coherence\result\cluster2.0\Cohere_freq')
a_vec2nii(brain0,'rank1');

cd('E:\GS Coherence\result\cluster2.0\Cohere_freq')
save Tmap_FDR T
a_vec2nii(T(1,:),'Tmap1-2.nii');
a_vec2nii(T(2,:),'Tmap1-3.nii');
a_vec2nii(T(3,:),'Tmap1-4.nii');
a_vec2nii(T(4,:),'Tmap2-3.nii');
a_vec2nii(T(5,:),'Tmap2-4.nii');
a_vec2nii(T(6,:),'Tmap3-4.nii');




