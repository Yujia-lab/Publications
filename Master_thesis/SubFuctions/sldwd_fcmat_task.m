function sldwd_fcmat_task(GS,GS_phase,ROISignals,wd_length,step_length,out_path,phase_amp)
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

imax = 360/step_length;

for isubj = 1:2:size(ROISignals,3)-1
    for iphase = 1:imax
        phase_tp1 = find_phase(GS_phase(:,isubj),iphase,wd_length,step_length);
        phase_tp2 = find_phase(GS_phase(:,isubj+1),iphase,wd_length,step_length);
        Comp_ROISignals = [ROISignals(phase_tp1,:,isubj);ROISignals(phase_tp2,:,isubj+1)];
        % compute GStopo
        fcmat(:,:,iphase) = corr(Comp_ROISignals,Comp_ROISignals,'type','pearson');
    end
    % all GS topo
    fcmat_all = corr(ROISignals(:,:,isubj),ROISignals(:,:,isubj),'type','pearson');
    %% save file
    cd([out_path, '\fcmat'])
    subj_id = num2str(isubj);
    save(['fcmat',subj_id],'fcmat')
    
    cd([out_path, '\fcmat_all'])
    save(['fcmat_all',subj_id],'fcmat_all')
    
    
end

    function phase_tp = find_phase(phase,iphase,wd_length,step_length)
        if iphase <= (360/step_length)-wd_length/step_length
            phase_tp = find(phase >= (-180+(iphase-1)*step_length) & phase <= (-180+wd_length+(iphase-1)*step_length));
        else
            phase_tp = find(phase >= (-180+(iphase-1)*step_length) | phase <= (-540+wd_length+(iphase-1)*step_length));
        end
    end

end