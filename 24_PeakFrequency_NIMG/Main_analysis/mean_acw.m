% calculate mean phase var across session
clear
clc


path = 'E:\Phase dynamics and INT\Results\acw_shuffled';
sub_file = dir(path);
sub_file(1:2) = [];

out_path = 'E:\Phase dynamics and INT\Results\mean_acw_shuffled';
mkdir(out_path)

for isub = 1:length(sub_file)
    data_file = dir(fullfile(path,sub_file(isub).name,'*.mat'));
    acw_all = [];
    
    for idata = 1:length(data_file)
        load(fullfile(path,sub_file(isub).name,data_file(idata).name));
        acw_all = cat(2,acw_all,acw_0'); 
    end

    mean_acw = mean(acw_all,2);
    name = sub_file(isub).name;
    id = name(1:6);
    cd(out_path)
    save(id,'mean_acw')
end
