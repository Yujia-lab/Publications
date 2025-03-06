%% phase variation for 3T
clc
clear


path = 'E:\Phase dynamics and INT\HCP_filt\3T';
ses_file = dir(fullfile(path,'*fMRI*'));

out_path = fullfile('E:\Phase dynamics and INT\Results','acw_shuffled');
mkdir(out_path)



for ises = 1:length(ses_file)
    ises

    path2 = fullfile(path,ses_file(ises).name);

    sub_file = dir(fullfile(path2,'*.mat'));


    for isub = 1:length(sub_file)
        isub
        load(fullfile(path2,sub_file(isub).name));

        for iROI = 1:size(roisignals_ff,2)
            data = squeeze(roisignals_ff(:,iROI));

            [acw_0(iROI), ~, acf, lags] = acw(data, 1/0.72);
        end

        sub_name = sub_file(isub).name;
        sub_name = sub_name(1:6);
        ses_name = ses_file(ises).name;
        ses_name = ses_name(end-7:end);
        name = [sub_name,'_',ses_name];
        mkdir(fullfile(out_path,sub_name))
        cd(fullfile(out_path,sub_name))
        save(name,'acw_0');

    end
end

