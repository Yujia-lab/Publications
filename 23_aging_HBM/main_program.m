clear
clc

%%  delete file 
%   Please use with caution, delete files can not be recovered!!!
%   usage:
%   start_path is the directory you need to delete
%   delete_name.txt is the file list you need to delete
clc;clear;
if exist('delete_name.txt','file')
    fid = fopen('delete_name.txt');
end
formatSpec = '%s';
delete_name = textscan(fid,formatSpec,'delimiter','\t');
fclose(fid);
start_path = uigetdir('start_path','please select a dir') ;
cd(start_path);
for i=1:numel(delete_name{1,1})
    if ~exist(delete_name{1,1}{i},'dir')==0
        rmdir(delete_name{1,1}{i}, 's');
    end
end

%%  计算GS
clear
clc
sub_path{1}='F:\lifespan\FunImgARWC_ALL_deconv';
sub_path{2}='F:\lifespan\FunImgARWC_ALL';
sub_path{3}='F:\lifespan\FunImgARWSC_Deconv\FunImgARWC';

mask_path='D:\program\matlab\Common\Mask\BN_Atlas_246_3mm.nii';    
result_path{1} = 'E:\GS Coherence\Deconv';
result_path{2} = 'E:\GS Coherence\Regression';
result_path{3} = 'E:\GS Coherence\Deconv_noreg';
for i = 1:3   
    i=3
    [ROISignals,GS] = a_GSROI(sub_path{i},mask_path);
    mkdir(result_path{i})
    cd(result_path{i});
    save GS GS
    save ROISignals ROISignals
end


%%   Coherence between ROIS and GS
clear
clc
sub_paths{1} = 'E:\GS Coherence\Deconv';
sub_paths{2} = 'E:\GS Coherence\Regression';
sub_paths{3} = 'E:\GS Coherence\Deconv_noreg';
%
TR=2;
nfft= 2^nextpow2(235*2);
window=hamming(32);
overlap=16;
fs=1/TR;
%
for i = 1:3
    GS_path = fullfile(sub_paths{i},'GS.mat');
    ROI_path = fullfile(sub_paths{i},'ROISignals.mat');
    load(GS_path);
    load(ROI_path);
    parfor isub = 1:492
        isub
        sub_ROISignals = ROISignals(:,:,isub);
        sub_GS = GS(:,isub);
        GStopo = zeros(257,246);
        FC = zeros(257,246,246);
        for iROI = 1:246            
            [c,f]= mscohere(sub_GS,sub_ROISignals(:,iROI),window,overlap,nfft,fs);
            ROI_Cohere(:,iROI,isub) = c;
        end
%         for iROI1 = 1:246   
%             for iROI2 = 1:246
%                 [c,f]= mscohere(sub_ROISignals(:,iROI1),sub_ROISignals(:,iROI2),window,overlap,nfft,fs);
%                 FC(:,iROI1,iROI2) = c;
%                 FC_Cohere{isub} = FC;
%             end
%         end
    end
    cd(sub_paths{i})
    Mean = mean(ROI_Cohere,3);
    save ROI_Cohere ROI_Cohere Mean
end

%% 相干单T检验
clear
clc
sub_paths{1} = 'E:\GS Coherence\Deconv';
sub_paths{2} = 'E:\GS Coherence\Regression';
sub_paths{3} = 'E:\GS Coherence\Non_regression';

for i = 1:3
    load(fullfile(sub_paths{i},'ROI_Cohere.mat'));
    for ifre = 1:257
        for iROI = 1:246
            Cohere = squeeze(ROI_Cohere(ifre,iROI,:));
            [h,p,ci,stats] = ttest(Cohere);
            T_Cohere(ifre,iROI) = stats.tstat;
            P_Cohere(ifre,iROI) = p;
        end
    end
end
            


%%   Correlation between the 'Cohere' and age
clear
clc
sub_paths{1} = 'E:\GS Coherence\Deconv';
sub_paths{2} = 'E:\GS Coherence\Regression';
sub_paths{3} = 'E:\GS Coherence\Deconv_noreg';

age_path='F:\lifespan\sub_information_492.xlsx';
info=xlsread(age_path);
N=length(info);
age=info(:,1);
gender=info(:,2); 
gender(gender==2) = 0;

for i = 1:3
    i
    load(fullfile(sub_paths{i},'ROI_Cohere.mat'))
    R = [];
    for iROI=1:246
        for iFreq=1:257
            X=[gender,ones(N,1)];%
            y=squeeze(ROI_Cohere(iFreq,iROI,:));
            [b,bint,Re,Rint,stats]=regress(y,X);   %
            [r,p]=corrcoef(Re,age);            
            R(iFreq,iROI) = r(2);
            P(iFreq,iROI) = p(2);
        end
    end
       
    Z = a_fishertrans(R);
    [m,n]=size(P);
    fdr=mafdr(reshape(P,m*n,1));
    label=find(fdr<0.05);
    R_FDR=zeros(m,n);
    R_FDR(label)=R(label);
    Z_FDR=zeros(m,n);
    Z_FDR(label)=Z(label);
    
    cd(sub_paths{i})
    save Age_GScorr R P R_FDR Z Z_FDR
end

%% correlation between GStopo and agetopo
clc
clear
data_path{1} = 'E:\GS Coherence\Deconv';
data_path{2} = 'E:\GS Coherence\Regression';
data_path{3} = 'E:\GS Coherence\Deconv_noreg';


for i = 1:3
    load(fullfile(data_path{i},'ROI_Cohere.mat'));
    load(fullfile(data_path{i},'Age_GScorr.mat'));   
    x = reshape(Mean,257*246,1);
    y = reshape(Z,257*246,1);
    [r,p] = corrcoef(x,y);
    R_all(i) = r(2);
    P_all(i) = p(2);
    
    %时间去分化
    for iROI = 1:246
        [r,p] = corrcoef(Mean(:,iROI),Z(:,iROI));
        R_temp(iROI) = r(2);
        P_temp(iROI) = p(2);
    end
    %空间去分化
    for ifre = 1:257
        [r,p] = corrcoef(Mean(ifre,:),Z(ifre,:));
        R_spatial(ifre) = r(2);
        P_spatial(ifre) = p(2);
    end

    % 分别计算时间去分化的正负R值
    R_tempfdr = a_multicorrect(R_temp,P_temp,'FDR');
    C_pos = Mean(:,R_tempfdr>0);
    C_neg = Mean(:,R_tempfdr<0);
    x2 = reshape(C_pos,[],1);
    x3 = reshape(C_neg,[],1);
    R_pos = Z(:,R_tempfdr>0);
    R_neg = Z(:,R_tempfdr<0);
    y2 = reshape(R_pos,[],1);
    y3 = reshape(R_neg,[],1);

    [r,p] = corrcoef(x2,y2);
    R_pos = r(2);
    P_pos = p(2);
    [r,p] = corrcoef(x3,y3);
    R_neg = r(2);
    P_neg = p(2);
    
    % 转ROI值为cifti文件方便可视化
    mask_path = 'D:\program\matlab\Common\Mask\BN_Atlas_freesurfer\fsaverage\fsaverage_LR164k\fsaverage.BN_Atlas.164k_fs_LR.dlabel.nii';
    cd(data_path{i})
    vec2cifti(R_tempfdr,mask_path,'temporal_dediff');

    % 分别计算空间去分化的正负R值
    R_spatialfdr = a_multicorrect(R_spatial,P_spatial,'FDR');
    C1 = mean(Mean(1:15,:),1);
    Z1 = mean(Z(1:15,:),1);
    C2 = mean(Mean(146:end,:),1);
    Z2 = mean(Z(146:end,:),1);
    [r,p] = corrcoef(C1,Z1);
    R1 = r(2);
    P1 = p(2);
    [r,p] = corrcoef(C2,Z2);
    R2 = r(2);
    P2 = p(2);
  

    cd(data_path{i})
    
    save spatial_dediff.mat R_temp P_temp R_spatial P_spatial R_tempfdr R_spatialfdr
    R_tempfdr(1:210)=0;
    a_vec2nii(R_tempfdr,'R_tempfdr.nii');
end

%% 做时空去分化的Z检验
clc
clear
data_path{1} = 'E:\GS Coherence\Deconv';
data_path{2} = 'E:\GS Coherence\Regression';
data_path{3} = 'E:\GS Coherence\Deconv_noreg';

load(fullfile(data_path{1},'spatial_dediff.mat'));
R_temp1 = R_temp;
R_spatial1 = R_spatial;
load(fullfile(data_path{2},'spatial_dediff.mat'));
R_temp2 = R_temp;
R_spatial2 = R_spatial;
load(fullfile(data_path{3},'spatial_dediff.mat'));
R_temp3 = R_temp;
R_spatial3 = R_spatial;

for i = 1:257
    [Zscore, Pvalue] = gretna_ztest_two_corrcoef(R_spatial2(i),R_spatial1(i), 257, 257, 'both');
    Z12(i) = Zscore;
    P12(i) = Pvalue;
    [Zscore, Pvalue] = gretna_ztest_two_corrcoef(R_spatial3(i),R_spatial1(i), 257, 257, 'both');
    Z13(i) = Zscore;
    P13(i) = Pvalue;
end
FDR_spaital12 = a_multicorrect(Z12,P12,'FDR');
FDR_spaital13 = a_multicorrect(Z13,P13,'FDR');

Z12 = [];
P12 = [];
Z13 = [];
P13 = [];
for i = 1:246
    [Zscore, Pvalue] = gretna_ztest_two_corrcoef(R_temp2(i), R_temp1(i), 257, 257, 'both');
    Z12(i) = Zscore;
    P12(i) = Pvalue;
    [Zscore, Pvalue] = gretna_ztest_two_corrcoef(R_temp3(i), R_temp1(i), 257, 257, 'both');
    Z13(i) = Zscore;
    P13(i) = Pvalue;
end
FDR_temp12 = a_multicorrect(Z12,P12,'FDR');
FDR_temp13 = a_multicorrect(Z13,P13,'FDR');
% 转ROI值为cifti文件方便可视化
mask_path = 'D:\program\matlab\Common\Mask\BN_Atlas_freesurfer\fsaverage\fsaverage_LR164k\fsaverage.BN_Atlas.164k_fs_LR.dlabel.nii';
cd(data_path{2})
vec2cifti(FDR_temp12,mask_path,'Z_test');
cd(data_path{3})
vec2cifti(FDR_temp13,mask_path,'Z_test');

% 可视化皮层下区域
cd(data_path{2})
FDR_temp12(1:210)=0;
a_vec2nii(FDR_temp12,'Z_sub.nii');

cd(data_path{3})
FDR_temp13(1:210)=0;
a_vec2nii(FDR_temp13,'Z_sub.nii');






