%% 3T
clear
clc

addpath('/home/yao/Downloads/matlab_toolbox/a_toolbox')
addpath('/home/yao/Downloads/matlab_toolbox/a_toolbox/Simulations-master')
addpath('/home/yao/Downloads/matlab_toolbox/cifti-matlab-master')

path = '/HCP/3T';
sub_file = dir(path);
sub_file(1:2) = [];
sub_file(end-1:end) = [];
tr = 0.72;
Fs = 1/tr;

template_path = '/home/yao/Downloads/matlab_toolbox/Glasser template/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
glasser_L = ciftiopen(template_path);
template_path = '/home/yao/Downloads/matlab_toolbox/Glasser template/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
glasser_R = ciftiopen(template_path);


for i = 1:length(sub_file)
    i

    cd(fullfile(path,sub_file(i).name));
    rest_folder = dir('*REST*fix');

    path2 = fullfile(path,sub_file(i).name,rest_folder.name,sub_file(i).name,'MNINonLinear','Results');
    data_file = dir(fullfile(path2,'rfMRI*'));

    if length(rest_folder) == 0
        continue
    end

    for i2 = 1:length(data_file)



        % check if the file exist
        file = dir(fullfile('/home/yao/Data/HCP/3T',sub_file(i).name,data_file(i2).name,'roisignals.mat'));
        if isempty(file) == 0
            continue
        end

        % load cifti and noise
        cifti_file = dir(fullfile(path2data_file(i2).name,'*MSMAll_hp2000_clean.dtseries.nii'));

        if length(cifti_file) == 0
            disp([sub_file(i).name,'no cifti'])
            continue          
        end

        surf = ciftiopen(fullfile(path2,data_file(i2).name,cifti_file.name));
        BOLD = surf.cdata;
        % extract roisignals for regressed data


        surf_L = BOLD(1:29696,:);
        surf_R = BOLD(29697:59412,:);


        label = unique(glasser_L.cdata);
        for iroi = 1:length(label)
            roi_pos = find(glasser_L.cdata == label(iroi));
            roi_voxel = surf_L(roi_pos,:);
            roisignals_L(iroi,:) = mean(roi_voxel,1);

            roi_pos = find(glasser_R.cdata == label(iroi));
            roi_voxel = surf_R(roi_pos,:);
            roisignals_R(iroi,:) = mean(roi_voxel,1);
        end

        roisignals = [roisignals_L;roisignals_R];


        mkdir(fullfile('/home/yao/Data/HCP/3T',sub_file(i).name,data_file(i2).name));
        cd(fullfile('/home/yao/Data/HCP/3T',sub_file(i).name,data_file(i2).name));
        save roisignals roisignals
        clear roisignals_L roisignals_R roisignals
    end
end

%% 7T
clear
clc

addpath('/home/yao/Downloads/matlab_toolbox/a_toolbox')
addpath('/home/yao/Downloads/matlab_toolbox/a_toolbox/Simulations-master')
addpath('/home/yao/Downloads/matlab_toolbox/cifti-matlab-master')

path = '/HCP/7T';
sub_file = dir(path);
sub_file(1:2) = [];
sub_file(end-1:end) = [];
tr = 0.72;
Fs = 1/tr;

template_path = '/home/yao/Downloads/matlab_toolbox/Glasser template/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
glasser_L = ciftiopen(template_path);
template_path = '/home/yao/Downloads/matlab_toolbox/Glasser template/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
glasser_R = ciftiopen(template_path);


for i = 1:length(sub_file)
    i

    cd(fullfile(path,sub_file(i).name));
    rest_folder = dir('*REST_2mm_fix');
    task_folder = dir('*MOVIE_2mm_fix');
    folder = [rest_folder;task_folder];

    if length(rest_folder) == 0
        continue
    end
    
    for i2 = 1:length(folder)
    
        path2 = fullfile(path,sub_file(i).name,folder(i2).name,sub_file(i).name,'MNINonLinear','Results');
        data_file = dir(fullfile(path2,'*fMRI*'));

        for i3 = 1:length(data_file)
         
        % load cifti and noise
        cifti_file = dir(fullfile(path2,data_file(i3).name,'*MSMAll_hp2000_clean.dtseries.nii'));
        
        if length(cifti_file) == 0
            disp([sub_file(i).name,'no cifti'])
            continue          
        end
        
        surf = ciftiopen(fullfile(path2,data_file(i3).name,cifti_file.name));
        
        
        BOLD = surf.cdata;
        
                   
        
        % extract roisignals for regressed data
        
        
        surf_L = BOLD(1:29696,:);
        surf_R = BOLD(29697:59412,:);


        label = unique(glasser_L.cdata);
        for iroi = 1:length(label)
            roi_pos = find(glasser_L.cdata == label(iroi));
            roi_voxel = surf_L(roi_pos,:);
            roisignals_L(iroi,:) = mean(roi_voxel,1);
            
            roi_pos = find(glasser_R.cdata == label(iroi));
            roi_voxel = surf_R(roi_pos,:);
            roisignals_R(iroi,:) = mean(roi_voxel,1);
        end
        
        roisignals = [roisignals_L;roisignals_R];
        
        
        if i2 == 1  
            mkdir(fullfile('/home/yao/Data/HCP/7T_rest',sub_file(i).name,data_file(i3).name));
            cd(fullfile('/home/yao/Data/HCP/7T_rest',sub_file(i).name,data_file(i3).name));
            save roisignals roisignals
            clear roisignals_L roisignals_R roisignals
        else
            mkdir(fullfile('/home/yao/Data/HCP/7T_task',sub_file(i).name,data_file(i3).name));
            cd(fullfile('/home/yao/Data/HCP/7T_task',sub_file(i).name,data_file(i3).name));
            save roisignals roisignals
            clear roisignals_L roisignals_R roisignals
        end
    end
    end
end

