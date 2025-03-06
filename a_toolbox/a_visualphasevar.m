function a_visualphasevar(mean_phase_var,path,name)

template_path1 = 'D:\Code\Common\Mask\Glasser_Template\Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
template_path2 = 'D:\Code\Common\Mask\Glasser_Template\Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
ciftipath = 'D:\Code\Common\Mask\BN_Atlas_freesurfer\fsaverage\fsaverage_LR32k\fsaverage.BN_Atlas.32k_fs_LR.dlabel.nii';
glasser_L = ciftiopen(template_path1);
glasser_R = ciftiopen(template_path2);
cifti = ciftiopen(ciftipath);

roi_data = squeeze(mean_phase_var);
for i = 1:length(roi_data)
    if i <= 180
        glasser_L.cdata(glasser_L.cdata == i) = roi_data(i);
    else
        glasser_R.cdata(glasser_R.cdata == i-180) = roi_data(i);
    end
end
cifti.cdata = [glasser_L.cdata;glasser_R.cdata];

mkdir(path)

cd(path)

cifti.diminfo{1,2}.type = 'scalars';
ciftisave(cifti,[name,'.dscalar.nii']);
