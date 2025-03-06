function [roisignals_ff,roisignals_filt] = a_butterfilt_revision(roisignals,Fs)

% filtering

low_cutoff = 0.01;
high_cutoff = 0.2;
bpFilt_ff = designfilt('bandpassiir', 'FilterOrder', 4,...
    'HalfPowerFrequency1', low_cutoff, 'HalfPowerFrequency2', high_cutoff,...
    'SampleRate', Fs, 'DesignMethod', 'butter');


low_cutoff = 0.09;
high_cutoff = 0.11;
bpFilt_fre1 = designfilt('bandpassiir', 'FilterOrder', 4,...
    'HalfPowerFrequency1', low_cutoff, 'HalfPowerFrequency2', high_cutoff,...
    'SampleRate', Fs, 'DesignMethod', 'butter');

low_cutoff = 0.11;
high_cutoff = 0.13;
bpFilt_fre2 = designfilt('bandpassiir', 'FilterOrder', 4,...
    'HalfPowerFrequency1', low_cutoff, 'HalfPowerFrequency2', high_cutoff,...
    'SampleRate', Fs, 'DesignMethod', 'butter');

low_cutoff = 0.13;
high_cutoff = 0.15;
bpFilt_fre3 = designfilt('bandpassiir', 'FilterOrder', 4,...
    'HalfPowerFrequency1', low_cutoff, 'HalfPowerFrequency2', high_cutoff,...
    'SampleRate', Fs, 'DesignMethod', 'butter');

low_cutoff = 0.15;
high_cutoff = 0.17;
bpFilt_fre4 = designfilt('bandpassiir', 'FilterOrder', 4,...
    'HalfPowerFrequency1', low_cutoff, 'HalfPowerFrequency2', high_cutoff,...
    'SampleRate', Fs, 'DesignMethod', 'butter');

low_cutoff = 0.17;
high_cutoff = 0.19;
bpFilt_fre5 = designfilt('bandpassiir', 'FilterOrder', 4,...
    'HalfPowerFrequency1', low_cutoff, 'HalfPowerFrequency2', high_cutoff,...
    'SampleRate', Fs, 'DesignMethod', 'butter');

parfor iroi = 1:size(roisignals,2)
    roidata = double(roisignals(:,iroi));
    roisignals_ff(:,iroi) = filtfilt(bpFilt_ff,roidata);

    for iband = 1:5

        if iband == 1
            bpFilt = bpFilt_fre1;
        elseif iband == 2
            bpFilt = bpFilt_fre2;
        elseif iband == 3
            bpFilt = bpFilt_fre3;
        elseif iband == 4
            bpFilt = bpFilt_fre4;
        elseif iband == 5
            bpFilt = bpFilt_fre5;
        end

        roisignals_filt(:,iroi,iband) = filtfilt(bpFilt,roisignals_ff(:,iroi));

    end
end