function [roisignals,GS] = a_GSROI(sub_path,mask_path)

% based on dpabi

sub_file = dir(sub_path);
sub_file(1:2) = [];
maskdata = y_Read(mask_path);
cd(sub_path);

parfor isub = 1:length(sub_file)
    isub
    
    
    % 读取被试的nii文件
    if sub_file(isub).isdir == 1
        brain_path = fullfile(sub_path, sub_file(isub).name);
        cd(brain_path);
        brain = dir('*.nii');
        [brain_data,header] = y_Read(brain.name);
    end
    
    if sub_file(isub).isdir == 0
        [brain_data,header] = y_Read(sub_file(isub).name); % read NII 
    end
    % 计算roisignal和GSsignal
    roisignals(:,:,isub) = a_roisignals(brain_data,maskdata);
    GS(:,isub) = a_globalsignal(brain_data,maskdata); % calculate GS  
end
