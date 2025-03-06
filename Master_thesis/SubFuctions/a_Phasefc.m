function FC = a_Phasefc(sub_path,mask_path,out_path,GS)

sub_file = dir(sub_path);
sub_file(1:2) = [];

maskdata=y_Read(mask_path);
dim = size(maskdata);
maskdata1=reshape(maskdata,dim(1)*dim(2)*dim(3), 1);
label0=find(maskdata==0);

for m = 1:4:length(sub_file)
    phase1 = a_Phase(GS(m,2:end));
    sub_GS1 = GS(m,2:end);
    sub_GS1(:,1:10) = [];
    sub_GS1(:,(end-9):end) = [];
    
    phase2 = a_Phase(GS(m+1,2:end));
    sub_GS2 = GS(m+1,2:end);
    sub_GS2(:,1:10) = [];
    sub_GS2(:,(end-9):end) = [];
    
    phase3 = a_Phase(GS(m+2,2:end));
    sub_GS3 = GS(m+2,2:end);
    sub_GS3(:,1:10) = [];
    sub_GS3(:,(end-9):end) = [];
    
    phase4 = a_Phase(GS(m+3,2:end));
    sub_GS4 = GS(m+3,2:end);
    sub_GS4(:,1:10) = [];
    sub_GS4(:,(end-9):end) = [];
    
    brain_path = fullfile(sub_path, sub_file(m).name);
    cd(brain_path);
    brain = dir('*.nii');
    [headdata1,header]=y_Read(brain.name);% 读取nii
    
    brain_path = fullfile(sub_path, sub_file(m+1).name);
    cd(brain_path);
    brain = dir('*.nii');
    [headdata2,header]=y_Read(brain.name);% 读取nii
    
    brain_path = fullfile(sub_path, sub_file(m+2).name);
    cd(brain_path);
    brain = dir('*.nii');
    [headdata3,header]=y_Read(brain.name);% 读取nii
    
    brain_path = fullfile(sub_path, sub_file(m+3).name);
    cd(brain_path);
    brain = dir('*.nii');
    [headdata4,header]=y_Read(brain.name);% 读取nii

    
    headdata1=reshape(headdata1,dim(1)*dim(2)*dim(3), length(headdata1));
    headdata1(:,(end-9):end) = [];
    headdata1(:,(1:10)) = [];
    
    headdata2=reshape(headdata2,dim(1)*dim(2)*dim(3), length(headdata2));
    headdata2(:,(end-9):end) = [];
    headdata2(:,(1:10)) = [];
    
    headdata3=reshape(headdata3,dim(1)*dim(2)*dim(3), length(headdata3));
    headdata3(:,(end-9):end) = [];
    headdata3(:,(1:10)) = [];
    
    headdata4=reshape(headdata4,dim(1)*dim(2)*dim(3), length(headdata4));
    headdata4(:,(end-9):end) = [];
    headdata4(:,(1:10)) = [];
    
    asde1 = find((phase1 >= -45 & phase1 <= 45)|(phase1 >= -180 & phase1 <= -135)|(phase1 >= 135 & phase1 <= 180));
    pt1 = find((phase1 >= 45 & phase1 <= 135)|(phase1 >= -135 & phase1 <= -45));
    
    asde2 = find((phase2 >= -45 & phase2 <= 45)|(phase2 >= -180 & phase2 <= -135)|(phase2 >= 135 & phase2 <= 180));
    pt2 = find((phase2 >= 45 & phase2 <= 135)|(phase2 >= -135 & phase2 <= -45));
    
    asde3 = find((phase3 >= -45 & phase3 <= 45)|(phase3 >= -180 & phase3 <= -135)|(phase3 >= 135 & phase3 <= 180));
    pt3 = find((phase3 >= 45 & phase3 <= 135)|(phase3 >= -135 & phase3 <= -45));
    
    asde4 = find((phase4 >= -45 & phase4 <= 45)|(phase4 >= -180 & phase4 <= -135)|(phase4 >= 135 & phase4 <= 180));
    pt4 = find((phase4 >= 45 & phase4 <= 135)|(phase4 >= -135 & phase4 <= -45));
    
    for iw = 1:246
        ROI_pos = find(maskdata1 == iw);
        ROI1_asde = headdata1(ROI_pos,asde1);
        ROI2_asde = headdata2(ROI_pos,asde2);
        ROI3_asde = headdata3(ROI_pos,asde3);
        ROI4_asde = headdata4(ROI_pos,asde4);
        ROI_asde = [ROI1_asde,ROI2_asde,ROI3_asde,ROI4_asde];
        
        ROI1_pt = headdata1(ROI_pos,pt1);
        ROI2_pt = headdata2(ROI_pos,pt2);
        ROI3_pt = headdata3(ROI_pos,pt3);
        ROI4_pt = headdata4(ROI_pos,pt4);
        ROI_pt = [ROI1_pt,ROI2_pt,ROI3_pt,ROI4_pt];
        
        tt = 1;
        for ivoxel1 = 1:length(ROI_pos)-1   
            ivoxel1
            for ivoxel2 = ivoxel1+1:length(ROI_pos)
                [r1,p1] = corrcoef(ROI_asde(ivoxel1,:),ROI_asde(ivoxel2,:));
                R1(tt) = r1(2);
                [r2,p2] = corrcoef(ROI_pt(ivoxel1,:),ROI_pt(ivoxel2,:));
                R2(tt) = r2(2);
                tt = tt+1;
            end
        end
        wFC_asde(m,iw) = mean(R1(:)); 
        wFC_pt(m,iw) = mean(R2(:));  
        RS_asde(iw,:) = mean(ROI_asde,1);
        RS_pt(iw,:) = mean(ROI_pt,1);
    end
    
    for iROI1 = 1:246
        iROI1
        for iROI2 = 1:246
            [r1,p1] = corrcoef(RS_asde(iROI1,:),RS_asde(iROI2,:));
            [r2,p2] = corrcoef(RS_pt(iROI1,:),RS_pt(iROI2,:));
            bFC_asde(iROI1,iROI2) = r1(2);
            bP_asde(iROI1,iROI2) = p1(2);
            bFC_pt(iROI1,iROI2) = r2(2);
            bP_pt(iROI1,iROI2) = p2(2);
        end
    end
    cd(out_path)
    mkdir([out_path,'/FC'])
    cd([out_path,'/FC'])
    name = [sub_file(m).name,'bFC.mat'];
    save(name,'bFC_asde','bP_asde','bFC_pt','bP_pt');
end
cd([out_path,'/FC'])
save wFC.mat wFC_asde wFC_pt