function vec2cifti(vector,ciftipath,name)

% convert vector to cifti file for brain map visulization

cifti = ciftiopen(ciftipath);
cifti.cdata(mod(cifti.cdata,2)==0 & cifti.cdata ~= 0) = cifti.cdata(mod(cifti.cdata,2)==0 & cifti.cdata ~= 0)-210;
for i = 1:210
    cifti.cdata(cifti.cdata == i) = vector(i);
end

cifti.diminfo{1,2}.type = 'scalars';
ciftisave(cifti,[name,'.dscalar.nii']);