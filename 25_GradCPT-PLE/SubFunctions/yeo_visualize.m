function yeo_visualize(roidata,path,name)

template_path = 'D:\Code\CBIG-master\stable_projects\brain_parcellation\Yeo2011_fcMRI_clustering\1000subjects_reference\Yeo_JNeurophysiol11_SplitLabels\fs_LR32k\Yeo2011_7Networks_N1000.dscalar.nii';
yeo_template = ciftiopen(template_path);

for i = 1:length(roidata)
        yeo_template.cdata(yeo_template.cdata == i) = roidata(i);
end

mkdir(path)
cd(path)
yeo_template.diminfo{1,2}.type = 'scalars';
ciftisave(yeo_template,[name,'.dscalar.nii']);