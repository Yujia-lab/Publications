function a_vec2nii(vector,name)

[~,header] = y_Read('/Users/yujia/Data/Data/Previous/lifespan_GS/3/0.08-0.1.nii');
mask_path = '/Users/yujia/Data/Code/Common/Mask/BN_Atlas_freesurfer/BN_Atlas_246_3mm.nii';
mask_data = y_Read(mask_path);

ROI_label = unique(mask_data);
iROI = 1:length(ROI_label)-1;
for i = iROI
    ROI_pos = find(mask_data == i);
    mask_data(ROI_pos) = vector(i);
end
y_Write(mask_data,header,name);