function [roisignals,GS] = a_GSROI(sub_path,mask_path)
%%% sub_path�Ǳ���·��
%%% mask_path��mask��nii�ļ���·��


sub_file = dir(sub_path);
sub_file(1:2) = [];

maskdata = y_Read(mask_path);
dim = size(maskdata);
label_0 = find(maskdata==0); % find mask label=0
cd(sub_path);

parfor isub = 1:length(sub_file)
    % ��ȡ���Ե�nii�ļ�
    if sub_file(isub).isdir == 1
        brain_path = fullfile(sub_path, sub_file(isub).name);
        cd(brain_path);
        brain = dir('*.nii');
        [brain_data,header] = y_Read(brain(1).name);
    end
    
    if sub_file(isub).isdir == 0
        [brain_data,header] = y_Read(sub_file(isub).name); % read NII 
    end
    % ����roisignal��GSsignal
    roisignals(:,:,isub) = a_roisignals(brain_data,maskdata);
    GS(:,isub) = a_globalsignal(brain_data,maskdata); % calculate GS  
end
