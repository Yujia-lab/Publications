clear
clc

subj_path = 'E:\GradCPT\BIDS_derivative';

subj_file = dir(fullfile(subj_path,'sub*'));

wbcmdPath = 'D:\Applications\workbench\bin_windows64\wb_command';

load('D:\Code\CBIG-master\stable_projects\brain_parcellation\Kong2019_MSHBM\lib\fs_LR_32k_medial_mask.mat');

template_path = 'D:\Code\CBIG-master\stable_projects\brain_parcellation\Yeo2011_fcMRI_clustering/1000subjects_reference\Yeo_JNeurophysiol11_SplitLabels\fs_LR32k\Yeo2011_7Networks_N1000.dlabel.nii';
yeo_network = ciftiopen(template_path);
label = 1:7;
surf_pos = yeo_network.cdata;

for isubj = 1:length(subj_file)
    isubj
    ses_file = dir(fullfile(subj_path,subj_file(isubj).name,'ses*'));
    for ises = 1:length(ses_file)
        cifti_file = dir(fullfile(subj_path,subj_file(isubj).name,ses_file(ises).name,'func','*.dtseries.nii'));

        cifti_path = fullfile(subj_path,subj_file(isubj).name,ses_file(ises).name,'func',cifti_file.name);
        ciftiData = ciftiopen(cifti_path,wbcmdPath);


        signals = ciftiData.cdata;
        signals = signals(1:59412,:);

        confound_file = dir(fullfile(subj_path,subj_file(isubj).name,ses_file(ises).name,'func','*.tsv'));
        confound_path = fullfile(subj_path,subj_file(isubj).name,ses_file(ises).name,'func',confound_file.name);


        confound = readtable(confound_path, "FileType","text",'Delimiter', '\t');


        csf = confound.csf;
        wm = confound.white_matter;
        FD = confound.framewise_displacement;
        hm1 = confound.trans_x;
        hm2 = confound.trans_y;
        hm3 = confound.trans_z;
        hm4 = confound.rot_x;
        hm5 = confound.rot_y;
        hm6 = confound.rot_z;


        x_csf = [csf, ones(size(hm1,1),1)];
        x_wm = [wm, ones(size(hm1,1),1)];
        x_FD = [FD, ones(size(hm1,1),1)];
        x_hm = [hm1,hm2,hm3,hm4,hm5,hm6, ones(size(hm1,1),1)];


        parfor ivoxel = 1:size(signals,1);

            [b,r,SSE,SSR] = y_regress_ss(signals(ivoxel,2:end)',x_csf(2:end,:));
            [b,r2,SSE,SSR] = y_regress_ss(r,x_wm(2:end,:));
            [b,r3,SSE,SSR] = y_regress_ss(r2,x_FD(2:end,:));
            [b,r4,SSE,SSR] = y_regress_ss(r3,x_hm(2:end,:));
            out_signals(ivoxel,:) = r4;
        end



        fMRI_data = zeros(64984,size(out_signals,2));
        fMRI_data(medial_mask,:) = out_signals;

        roisignals = zeros(8,size(out_signals,2));

        roisignals(1,:) = mean(out_signals,1);

        for iroi = 1:length(label)
            roi_pos = find(surf_pos == label(iroi));
            roi_voxel = fMRI_data(roi_pos,:);
            roisignals(iroi+1,:) = mean(roi_voxel,1);
        end

        path3 = fullfile('E:\GradCPT\network_signals',ses_file(ises).name);
        mkdir(path3);
        cd(path3)
        name = subj_file(isubj).name;
        save(name,'roisignals');

        clear roi_pos roi_voxel out_signals r r2 r3 r4

    end
end