% calculate mean phase var across session
clear
clc

dataset = '7T_task'
path = ['E:\Phase dynamics and INT\Results\phase_mean_highfrequency_',dataset];
sub_file = dir(path);
sub_file(1:2) = [];

out_path = ['E:\Phase dynamics and INT\Results\mean_phase_highfrequency_',dataset];
mkdir(out_path)

for isub = 1:length(sub_file)
    isub
    data_file = dir(fullfile(path,sub_file(isub).name,'*.mat'));
    phase_mean_dy_all = [];
    phase_mean_stat_all = [];
    phase_mean_dy_shuffled_all = [];
    phase_mean_stat_shuffled_all = [];

    for idata = 1:length(data_file)
        load(fullfile(path,sub_file(isub).name,data_file(idata).name));

        phase_mean_dy_all = cat(4,phase_mean_dy_all,phase_mean_dy);       
        phase_mean_stat_all = cat(3,phase_mean_stat_all,phase_mean_stat);

        % phase_mean_dy_shuffled_all = cat(4,phase_mean_dy_shuffled_all,phase_mean_dy_shuffled);       
        % phase_mean_stat_shuffled_all = cat(3,phase_mean_stat_shuffled_all,phase_mean_stat_shulffed);
    end

    mean_phase_mean_dy = mean(phase_mean_dy_all,4);
    mean_phase_mean_stat = mean(phase_mean_stat_all,3);

    % mean_phase_mean_dy_shuffled = mean(phase_mean_dy_shuffled_all,4);
    % mean_phase_mean_stat_shuffled = mean(phase_mean_stat_shuffled_all,3);


    name = sub_file(isub).name;
    id = name(1:6);
    cd(out_path)
    save(id,'mean_phase_mean_stat','mean_phase_mean_dy')
end







