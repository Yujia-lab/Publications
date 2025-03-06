%% phase variation for 3T
clc
clear
dataset{1} = '3T';
dataset{2} = '7T_rest';
dataset{3} = '7T_task';

step_length = 1;
wd_length = 30;
imax = 360/step_length;
for iset = 1:1

    path = fullfile('E:\Phase dynamics and INT\HCP_filt',dataset{iset});
    ses_file = dir(fullfile(path,'*fMRI*'));

    out_path = fullfile('E:\Phase dynamics and INT\Results',['fp_',dataset{iset}]);
    mkdir(out_path)

    for ises = 1:length(ses_file)
        ises

        path2 = fullfile(path,ses_file(ises).name);

        sub_file = dir(fullfile(path2,'*.mat'));


        for isub = 1:length(sub_file)
            isub
            load(fullfile(path2,sub_file(isub).name));

            parfor iband = 1:size(roisignals_filt,3)
                for iROI = 1:360

                    data = squeeze(roisignals_filt(:,iROI,iband));
                    data_hil = hilbert(data);
                    realPart = real(data_hil);
                    imgPart = imag(data_hil);

                    fp = sqrt(realPart.^2+imgPart.^2);
                    fp(1:10) = [];
                    fp(end-9:end) = [];

                    fp_mean(iROI,iband) = mean(fp);
                    data_phase2 = a_Phase(data);

                    for iphase = 1:imax

                        if iphase <= (360/step_length)-wd_length/step_length
                            phase_tp = find(data_phase2 >= (-180+(iphase-1)*step_length) & data_phase2 <= (-180+wd_length+(iphase-1)*step_length));
                        else
                            phase_tp = find(data_phase2 >= (-180+(iphase-1)*step_length) | data_phase2 <= (-540+wd_length+(iphase-1)*step_length));
                        end


                        fp_dy(iROI,iband,iphase) = mean(fp(phase_tp));

                    end

           




                end
            end

            fp_dy = cat(3,fp_dy(:,:,end-14:end),fp_dy(:,:,1:end-15));




            sub_name = sub_file(isub).name;
            sub_name = sub_name(1:6);
            ses_name = ses_file(ises).name;
            ses_name = ses_name(end-7:end);
            name = [sub_name,'_',ses_name];
            mkdir(fullfile(out_path,sub_name))
            cd(fullfile(out_path,sub_name))
            save(name,'fp_mean','fp_dy');


        end
    end
end
