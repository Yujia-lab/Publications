function sldwd_fcmat(GS,GS_phase,ROISignals,wd_length,step_length,out_path,phase_amp)
GS_phase(1:10,:) = [];
GS_phase(end-9:end,:) = [];
ROISignals(1:10,:,:) = [];
ROISignals(end-9:end,:,:) = [];
GS(1:10,:) = [];
GS(end-9:end,:) = [];

mkdir([out_path, '\fcmat_all'])
mkdir([out_path, '\fcmat'])
mkdir([out_path, '\fcmat_all_amp'])
mkdir([out_path, '\fcmat_amp'])
if phase_amp == 0
    imax = 360/step_length;
elseif phase_amp == 1
    imax = floor((length(GS)-round(size(ROISignals,1)*wd_length/100))/(round(size(ROISignals,1)*step_length/100)))+1;
end
for isubj = 1:size(ROISignals,3)
    for iphase = 1:imax
        if phase_amp == 0
            phase_tp = find_phase(GS_phase(:,isubj),iphase,wd_length,step_length);
        elseif phase_amp == 1
            phase_tp = find_amp(GS(:,isubj),iphase,wd_length,step_length);
        end
        % compute GStopo
        fcmat(:,:,iphase) = corr(ROISignals(phase_tp,:,isubj),ROISignals(phase_tp,:,isubj),'type','pearson');
    end
    % all GS topo
    fcmat_all = corr(ROISignals(:,:,isubj),ROISignals(:,:,isubj),'type','pearson');
    %% save file
    if phase_amp == 0
        cd([out_path, '\fcmat'])
        subj_id = num2str(isubj);
        save(['fcmat',subj_id],'fcmat')
        
        cd([out_path, '\fcmat_all'])
        save(['fcmat_all',subj_id],'fcmat_all')
    elseif phase_amp == 1
        cd([out_path, '\fcmat_amp'])
        subj_id = num2str(isubj);
        save(['fcmat',subj_id],'fcmat')
        
        cd([out_path, '\fcmat_all_amp'])
        save(['fcmat_all',subj_id],'fcmat_all')
    end
    
end

    function phase_tp = find_phase(phase,iphase,wd_length,step_length)
        if iphase <= (360/step_length)-wd_length/step_length
            phase_tp = find(phase >= (-180+(iphase-1)*step_length) & phase <= (-180+wd_length+(iphase-1)*step_length));
        else
            phase_tp = find(phase >= (-180+(iphase-1)*step_length) | phase <= (-540+wd_length+(iphase-1)*step_length));
        end
    end

    function phase_tp = find_amp(GS,i,wd_length,step_length)
        [B,I] = sort(GS,'descend');
        N = round(length(I)*wd_length/100);
        N2 = round(length(I)*step_length/100);
        max = (1+(i-1)*N2);
        min = (1+(i-1)*N2)+N;
        phase_tp = I(max:min);
    end
end