function comp_phase_diff(phase_diff,out_path)
absvalue = abs(phase_diff);
for i = 1:246
    [value1,pos1] = max(absvalue(50:100,i));
    pos1 = pos1 + 49;
    maxvalue(i) = phase_diff(pos1,i);
    [value2,pos2] = max(absvalue(220:260,i));
    pos2 = pos2 + 219;
    all_pos1(i) = pos1;
end

mask_path = 'D:\program\matlab\Common\Mask\BN_Atlas_246_2mm.nii';

maskdata = y_Read(mask_path);
for iROI = 1:246
    ROI_pos = find(maskdata == iROI);
    maskdata(ROI_pos) = maxvalue(iROI);
end
cd(out_path)
y_Write(maskdata,header,'phase_diff');
