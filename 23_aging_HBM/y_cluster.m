clear
clc

%%   Import a header file
brain_path='D:\yangchengxiao\DATA_lifespan\lifespan\FunRawARWSC\sub31275\CovRegressed_4DVolume.nii';
[brain_data,header]=y_Read(brain_path);
%%  Import a mask file
mask_path='D:\yangchengxiao\matlab-aoyujia\Common\Mask\BN_Atlas_246_3mm.nii';
mask_data=y_Read(mask_path);
%%  Import data file--cohere
data_path1='D:\yangchengxiao\DATA_lifespan2.0\result\ROIS_cohere';
cd(data_path1);
load('GS_topo_Cohere.mat');
data=mean(ROI_Cohere,3);       %   Average by subject
%%   Import data file--corr
data_path2='D:\yangchengxiao\DATA_lifespan2.0\result\CORR';
cd(data_path2);
load('GSCORR.mat');
data=Z;               %   


%%  Binarization of the correlation value (Z) :
%   FDR>0 assigns a value of 1, <0 assigns a value of -1, and the others are 0
data_path='D:\yangchengxiao\DATA_lifespan2.0\result\CORR';
cd(data_path);
load('GSCORR.mat');

[m,n]=size(Z);
Z_Binary=zeros(m,n);
for i=1:m
    for j=1:n
        if FDR(i,j)>0
            Z_Binary(i,j)=1;
        elseif FDR(i,j)<0
                Z_Binary(i,j)=-1;
        else
            Z_Binary(i,j)=0;
        end      
    end
end

data=Z_Binary;
save GSCORR_Binary1.mat Z_Binary


%%  Binarization of the correlation value (FDR) :
%   90% threshold binarization
data_path='D:\yangchengxiao\DATA_lifespan2.0\result\CORR';
cd(data_path);
load('GSCORR_zscores.mat');

[m,n]=size(FDR);
Z_Binary=zeros(m,n);
P_sort=sort(P,2,'descend');
label_top=round(n*0.9);
PValue_label=P_sort(:,label_top);
for i=1:m
    for j=1:n
        if P(i,j)<=PValue_label(i)
            if Z(i,j)>0
               Z_Binary(i,j)=1;
            elseif Z(i,j)<0
               Z_Binary(i,j)=-1;
            else
               Z_Binary(i,j)=0;
            end
        else    
            Z_Binary(i,j)=0;
        end      
    end
end

data=Z_Binary;
save GSCORR_Binary2.mat Z_Binary

%%   Estimate the number of clusters
obj = evalclusters(data, 'kmeans','CalinskiHarabasz','klist',2:8);

% k_num=12;
% bestD = evaclusters_elbow(data,k_num);
%%   kmeans cluster
num_cluster=obj.OptimalK;
num_cluster=3;
opts = statset('MaxIter', 500, 'Display', 'off');
[IDX, C, SUMD, D] = kmeans(data, num_cluster,'Replicates',5,'options',opts);

%%   Clustering results of brain regions,Convert to NII file and save
name=strcat('cluster_region');      %   
vector=IDX;
result_path='D:\yangchengxiao\DATA_lifespan2.0\result\cluster2.0\CORR\region';
cd(result_path);
a_vec2nii(vector,mask_data,header,name)
save cluster_region.mat data IDX C

%%   Calculate the sequence of IDX, used to calculate the frequency band
sequence=zeros(1,num_cluster);          %   
AA=zeros(1,num_cluster);                %   
number=1;
a=1;
for i=1:length(IDX)-1
    if IDX(i)==IDX(i+1)
        number=number+1;
        AA(a)=number;
        sequence(a)=IDX(i);
    else  
        number=1;
        a=a+1; 
    end 
end
%%   Calculate the frequency range of the cluster
n_band=length(sequence);
bandrange=zeros(n_band,2);
band1=zeros(n_band,1);
band2=zeros(n_band,1);
for j=1:n_band
    n_freq_point=AA(j);
    band2(j)=band1(j)+0.25/500.*n_freq_point;
    k=j+1;
    if k<=n_band
    band1(k)=band2(j);
    else
    end
end
    bandrange(:,1)=band1;
    bandrange(:,2)=band2;

%%   Brain maps of each cluster category
[x,y]=size(data);
mean_data=zeros(n_band,y);  
AAA=[0,AA];
iii=1;
for ii=1:n_band
    c=sum(AAA(1:iii))+1;
    d=sum(AAA(1:iii+1));
    mean_data(ii,:)=mean(data(c:d,:),1);     %   The average value of each cluster band
    iii=iii+1;

    mask_data=y_Read(mask_path);
    ROI_label=unique(mask_data);
    n=length(ROI_label)-1; 
    for jj=1:n
        ROI_pos=find(mask_data==jj);
        mask_data(ROI_pos)=mean_data(ii,jj);
    end
    
name=strcat('CORR_4-',num2str(ii));      %   Create filename
result_path='D:\yangchengxiao\DATA_lifespan2.0\result\cluster2.0\CORR\freq';
cd(result_path);
y_Write(mask_data,header,name);
end
save cluster_frequency.mat mean_data IDX sequence AA bandrange     %   save

%%   Brain maps of each cluster category
%    Z values transfer to 'r' to 'P' and then FDR correction
[x,y]=size(data);
mean_data=zeros(n_band,y);  
AAA=[0,AA];
iii=1;
for ii=1:n_band
    c=sum(AAA(1:iii))+1;
    d=sum(AAA(1:iii+1));
    mean_data(ii,:)=mean(data(c:d,:),1);     %   The average value of each cluster band
    iii=iii+1;
end
%   Fisher transformation, evaluate 'r'
r_of_z= gretna_inv_fishertrans(mean_data);

rv=r_of_z;
nv=322;
PP=a_r2p(rv,nv);   %   'r' to P

%   FDR correction
[m,n]=size(PP);
fdr=mafdr(reshape(PP,m*n,1));
FDR_Z=zeros(m,n);
position=find(fdr<0.05);
FDR_Z(position)=mean_data(position);

for jj=1:n_band
    mask_data=y_Read(mask_path);
    ROI_label=unique(mask_data);
    n=length(ROI_label)-1; 
    for jjj=1:n
        ROI_pos=find(mask_data==jjj);
        mask_data(ROI_pos)=FDR_Z(jj,jjj);
    end
    
name=strcat('CORR_4-',num2str(jj));      %   Create filename
result_path='D:\yangchengxiao\DATA_lifespan2.0\result\cluster2.0\CORR\freq';
cd(result_path);
y_Write(mask_data,header,name);
end
save cluster_frequency.mat mean_data FDR_Z IDX sequence AA bandrange     %   save
