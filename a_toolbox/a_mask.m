function mask_brain = a_mask(brain_data, mask_data);


dim = size(mask_data);
label_0 = find(mask_data==0);
mask_brain=reshape(brain_data, dim(1) * dim(2) * dim(3), size(brain_data,4));
mask_brain(label_0,:)=[];
