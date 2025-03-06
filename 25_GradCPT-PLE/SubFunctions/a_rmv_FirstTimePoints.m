function a_rmv_FirstTimePoints(Indir,TimePoints)
% remove n first time points for BIDS data, replace the original data
% only support for nifti data

subj_folder = dir(fullfile(Indir,'sub*'));

for isubj = 1:length(subj_folder)
    img_file = dir(fullfile(Indir,subj_folder(isubj).name,'func','*.nii'));
    disp(['removing First Time Points for ',subj_folder(isubj).name,'...'])
    DirImg = fullfile(Indir,subj_folder(isubj).name,'func',img_file.name);
    [Data Header]=y_Read(DirImg);
    y_Write(Data(:,:,:,TimePoints+1:end),Header,DirImg);
end

