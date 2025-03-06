function sldwd_GStopo(GS,GS_phase,ROISignals,wd_length,step_length,out_path)

GS_phase(1:10,:) = [];
GS_phase(end-9:end,:) = [];
GS(1:10,:) = [];
GS(end-9:end,:) = [];
ROISignals(1:10,:,:) = [];
ROISignals(end-9:end,:,:) = []; %去前后十个时间点
mkdir([out_path, '\sub_topo'])
mkdir([out_path, '\sub_topo_all'])
mkdir([out_path, '\sub_topo_amp'])
mkdir([out_path, '\sub_topo_all_amp']) %创建目录

imax = 360/step_length;

for isubj = 1:size(GS,2)
    for iphase = 1:imax
        
        phase_tp = find_phase(GS_phase(:,isubj),iphase,wd_length,step_length);

        % compute GStopo
        GS_topo(:,iphase) = corr(ROISignals(phase_tp,:,isubj),GS(phase_tp,isubj),'type','pearson');
    end
    % all GS topo
    GS_topo_all = corr(ROISignals(:,:,isubj),GS(:,isubj),'type','pearson');
    %% save file
       

    cd([out_path, '\sub_topo'])
    subj_id = num2str(isubj);
    save(['GS_topo',subj_id],'GS_topo')


    cd([out_path, '\sub_topo_all'])
    save(['GS_topo_all',subj_id],'GS_topo_all')

end

    function phase_tp = find_phase(phase,iphase,wd_length,step_length)
        if iphase <= (360/step_length)-wd_length/step_length
            phase_tp = find(phase >= (-180+(iphase-1)*step_length) & phase <= (-180+wd_length+(iphase-1)*step_length));
        else
            phase_tp = find(phase >= (-180+(iphase-1)*step_length) | phase <= (-540+wd_length+(iphase-1)*step_length));
        end
    end

end