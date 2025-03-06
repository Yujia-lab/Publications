%% phase variation for 3T
clc
clear
dataset{1} = '3T';
dataset{2} = '7T_rest';
dataset{3} = '7T_task';


for iset = 1:1

    path = fullfile('E:\Phase dynamics and INT\HCP_filt',dataset{iset});
    ses_file = dir(fullfile(path,'*fMRI*'));

    out_path = fullfile('E:\Phase dynamics and INT\Results',['phase_mean_',dataset{iset}]);
    mkdir(out_path)

    for ises = 1:length(ses_file)
        ises

        path2 = fullfile(path,ses_file(ises).name);

        sub_file = dir(fullfile(path2,'*.mat'));


        for isub = 1:length(sub_file)
            isub
            load(fullfile(path2,sub_file(isub).name));

            step_length = 1;
            wd_length = 30;
            imax = 360/step_length;

            parfor iband = 1:size(roisignals_filt,3)
                for iROI = 1:360

                    data = squeeze(roisignals_filt(:,iROI,iband));
                    data_hil = hilbert(data);
                    data_phase = phase(data_hil);
                    data_phase(1:10) = [];
                    data_phase(end-9:end) = [];
                    data_phase2 = a_Phase(data);

                    data_diff = diff(data_phase);

                    filt_data_diff = medfilt1(data_diff,20,'omitnan','truncate');
                    phase_mean_stat(iROI,iband) = mean(filt_data_diff);
                    for iphase = 1:imax

                        if iphase <= (360/step_length)-wd_length/step_length
                            phase_tp = find(data_phase2 >= (-180+(iphase-1)*step_length) & data_phase2 <= (-180+wd_length+(iphase-1)*step_length))-1;
                        else
                            phase_tp = find(data_phase2 >= (-180+(iphase-1)*step_length) | data_phase2 <= (-540+wd_length+(iphase-1)*step_length))-1;
                        end

                        phase_tp(phase_tp==0) = [];
                        phase_mean_dy(iROI,iband,iphase) = mean(filt_data_diff(phase_tp));

                    end


                    data = squeeze(roisignals_filt(:,iROI,iband));
                    data = phaseran(data,1);
                    data_hil = hilbert(data);
                    data_phase = phase(data_hil);
                    data_phase(1:10) = [];
                    data_phase(end-9:end) = [];
                    data_phase2 = a_Phase(data);

                    data_diff = diff(data_phase);

                    filt_data_diff = medfilt1(data_diff,20,'omitnan','truncate');
                    phase_mean_stat_shulffed(iROI,iband) = mean(filt_data_diff);
                    for iphase = 1:imax

                        if iphase <= (360/step_length)-wd_length/step_length
                            phase_tp = find(data_phase2 >= (-180+(iphase-1)*step_length) & data_phase2 <= (-180+wd_length+(iphase-1)*step_length))-1;
                        else
                            phase_tp = find(data_phase2 >= (-180+(iphase-1)*step_length) | data_phase2 <= (-540+wd_length+(iphase-1)*step_length))-1;
                        end

                        phase_tp(phase_tp==0) = [];
                        phase_mean_dy_shuffled(iROI,iband,iphase) = mean(filt_data_diff(phase_tp));

                    end
                end
            end


            phase_mean_dy = cat(3,phase_mean_dy(:,:,end-14:end),phase_mean_dy(:,:,1:end-15));
            phase_mean_dy_shuffled = cat(3,phase_mean_dy_shuffled(:,:,end-14:end),phase_mean_dy_shuffled(:,:,1:end-15));


            sub_name = sub_file(isub).name;
            sub_name = sub_name(1:6);
            ses_name = ses_file(ises).name;
            ses_name = ses_name(end-7:end);
            name = [sub_name,'_',ses_name];
            mkdir(fullfile(out_path,sub_name))
            cd(fullfile(out_path,sub_name))
            save(name,'phase_mean_dy','phase_mean_stat','phase_mean_stat_shulffed','phase_mean_dy_shuffled');
            clear phase_mean_dy_shuffled

        end
    end
end
