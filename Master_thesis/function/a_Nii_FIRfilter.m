function a_Nii_FIRfilter(sub_path, out_path, bands, Fs)
%%% Coding by Aoyujia 2020.7.20
%%% See manu for more instruction


sub_file = dir(sub_path);
sub_file(1:2) = [];

mkdir(out_path)
cd(out_path)


 for m = 1:length(bands)

    
    DEV = [0.05 0.01 0.05];
    A = [0 1 0];
    band = bands{m}
    F = [band(1)-0.01, band(1), band(2), band(2)+0.01];
    [N, Fo, Ao, W] = firpmord(F, A, DEV, Fs);
    if(mod(N, 2) == 1) 
        N = N + 1; % it has to be even -> b is odd -> delay is integer = N/2
    end
    psess.filter.b = firpm(N, Fo, Ao, W);
    psess.filter.N = N;
    
    out_path2 = fullfile(out_path, num2str(band));
    mkdir(out_path2);
    
     parfor i = 1:length(sub_file)
        if sub_file(i).isdir == 1
            brain_path = fullfile(sub_path, sub_file(i).name);
            cd(brain_path);
            brain = dir('*.nii');
            [brain_data,header]=y_Read(brain.name);% ∂¡»°nii
        end

        if sub_file(i).isdir == 0
            [brain_data,header]=y_Read(sub_file(i).name);
        end
        
        disp(sub_file(i).name)
        dim = size(brain_data);
        brain_data = reshape(brain_data, dim(1)*dim(2)*dim(3), dim(4));
        
        for j = 1:dim(1)*dim(2)*dim(3)
            if max(brain_data(j, :)) ~= min(brain_data(j, :))  
                ts = brain_data(j, :);
                tsnorm = ts - mean(ts);             % removing DC for optimal bandpass filtering
                T = length(tsnorm);
                tsout = conv(psess.filter.b, tsnorm);  % produces a T+N-1 length signal
                tsout = flipud(conv(psess.filter.b, flipud(tsout)));
                if(0) % delay compensation, normal filter [not implemented]
                    tsout = tsout((psess.filter.N/2 + 1):(end - psess.filter.N/2));
                else
                    tsout = tsout(length(psess.filter.b):end);
                    tsout(T+1:end) = [];
                end
                if(0) % option for filter transient removal
                    tsout = tsout((psess.filter.N + 1):end); % removing N filter transient and N/2 filter tail of factual data
                end
                brain_data(j, :) = tsout;
            end
        end
        
        cd(out_path2) 
        brain_data = reshape(brain_data, dim(1), dim(2), dim(3), dim(4));  
        y_Write(brain_data,header,[sub_file(i).name,'.nii'])
    end
 end

        