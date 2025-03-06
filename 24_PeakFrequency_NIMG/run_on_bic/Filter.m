%% For 3T
clc
clear
subj_path = 'E:\Data\HCP\HCP\3T';
subj_file = dir(subj_path);
subj_file(1:2) = [];

Fs = 1/0.72;

for i = 1:length(subj_file)
    try
    i
    data_file = dir(fullfile(subj_path,subj_file(i).name,'rfMRI*'));
    for i2 = 1:length(data_file)
                
        load(fullfile(subj_path,subj_file(i).name,data_file(i2).name,'roisignals.mat'));
               
        % filtering
        [roisignals_ff,roisignals_filt] = a_butterfilt(roisignals',Fs);
               
        mkdir(fullfile('E:\Data\HCP\HCP_filt\3T',data_file(i2).name))
        cd(fullfile('E:\Data\HCP\HCP_filt\3T',data_file(i2).name))
        name = subj_file(i).name;
        save(name,'roisignals_ff','roisignals_filt')
        clear roisignals_ff  roisignals_filt
    end
    end
end


%% For 7T rest
clc
clear

subj_path = 'E:\Data\HCP\HCP\7T_rest';
subj_file = dir(subj_path);
subj_file(1:2) = [];

Fs = 1/0.72;


for i = 1:length(subj_file)
    try
    i
    data_file = dir(fullfile(subj_path,subj_file(i).name));
    data_file(1:2) = [];
    for i2 = 1:length(data_file)
        
                
        load(fullfile(subj_path,subj_file(i).name,data_file(i2).name,'roisignals.mat'));
               
        % filtering
        [roisignals_ff,roisignals_filt] = a_butterfilt(roisignals',Fs);
        
        
        mkdir(fullfile('E:\Data\HCP\HCP_filt\7T_rest',data_file(i2).name))
        cd(fullfile('E:\Data\HCP\HCP_filt\7T_rest',data_file(i2).name))
        name = subj_file(i).name;
        save(name,'roisignals_ff','roisignals_filt')
        clear roisignals_ff  roisignals_filt
        end
    end
end

%% For 7T task
clc
clear

subj_path = 'E:\Data\HCP\HCP\7T_task';
subj_file = dir(subj_path);
subj_file(1:2) = [];

Fs = 1/0.72;




for i = 1:length(subj_file)
    try
    i
    data_file = dir(fullfile(subj_path,subj_file(i).name));
    data_file(1:2) = [];
    for i2 = 1:length(data_file)
        
        
        load(fullfile(subj_path,subj_file(i).name,data_file(i2).name,'roisignals.mat'));
               
        % filtering
        [roisignals_ff,roisignals_filt] = a_butterfilt(roisignals',Fs);
        
        
        mkdir(fullfile('E:\Data\HCP\HCP_filt\7T_task',data_file(i2).name))
        cd(fullfile('E:\Data\HCP\HCP_filt\7T_task',data_file(i2).name))
        name = subj_file(i).name;
        save(name,'roisignals_ff','roisignals_filt')
        clear roisignals_ff  roisignals_filt
    end
    end
end

