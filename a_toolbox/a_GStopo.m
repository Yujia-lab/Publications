function a_GStopo(sub_path, out_path, mask_path)
%%% Coding by Aoyujia 2020.7.20
%%% See manu for more instruction

%% read dir
sub_file = dir(sub_path);
sub_file(1:2) = [];

maskdata=y_Read(mask_path);
dim = size(maskdata);
label1=find(maskdata~=0);
label0=find(maskdata==0);

cd(sub_path);


%% compute GS topography

for m=1:length(sub_file)
    m
    % read NII data
    if  sub_file(m).isdir == 1
        brain_path = fullfile(sub_path, sub_file(m).name);
        cd(brain_path);
        brain = dir('*.nii');
        [headdata,header]=y_Read(brain.name);% ∂¡»°nii
    end
    
    if sub_file(m).isdir == 0
        [headdata,header]=y_Read(sub_file(m).name);
    end
    
    headdata1=reshape(headdata,dim(1)*dim(2)*dim(3), length(headdata));
    headdata1(label0,:)=[];
    R = [];
    Z = [];
    GS = mean(headdata1,1);
    
    
    % compute correlation
    for i=1:length(headdata1)
        r = corrcoef(headdata1(i,:), GS);
        R(1,i)=r(2);
        Z(1,i) = 0.5*log((1+r(2))/(1-r(2)));
    end
    
    all_R(m, :) = R;
    all_Z(m, :) = Z;
end

%% save file
cd(out_path);

% all GS topo save in mat
save('all_topo.mat','all_R');
save(['all_topo_Z.mat'],'all_Z');
mkdir([out_path, '\subj_topo']);
cd([out_path, '\subj_topo']);

% GS topo for each subj
for n = 1:length(sub_file)
    a_mat2nii(all_Z(n,:), maskdata,header,[sub_file(n).name]); 
end

% mean GS topo
meandata1 = mean(all_R);
meandata2 = mean(all_Z);
cd(out_path)
a_mat2nii(meandata1, maskdata,header,'R_topo');
a_mat2nii(meandata2, maskdata,header,'Z_topo');