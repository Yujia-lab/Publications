function [PLE, psd, freq] = ple(ts, fs, freqrange, isplot)
%% Power-Law Exponent
% Neurophysiological measurements often have scale-free power spectrum. This
% power spectrum is characterized by the relation psd = 1 / freq^alpha.
% Here, alpha (> 0) measures how fast the power decreases with increasing
% frequency. An alpha close to 1 might indicate (self-organized) criticality: being
% close to phase transitions (think ice becoming water when temperature is
% increased from -1 to 1). 
% These two links are great resources:
% http://www.scholarpedia.org/article/Scale-free_neocortical_dynamics
% http://www.scholarpedia.org/article/1/f_noise
% For more thorough treatments, see the literature of Bruce West.
%
% It should be noted that linking 1/f noise to self-organization or making
% any inference whatsoever is controversial. Many other distributions 
% can look like power-laws locally. For further reference, see 
% Stumpf and Porter, 2012 (https://doi.org/10.1126/science.1216142) and 
% Clauset, Shalizi and Newman, 2009 (http://arxiv.org/abs/0706.1062).

% Input:
%   ts: time series (1D vector)
%   fs: sampling frequency (scalar)
%   freqrange: frequency range you are interested in ([fmin fmax]). fmin >
%   0 (since log(0) = inf).
%   isplot: plot to see the PSD? (logical)
% Output: 
%   PLE: power-law exponent
%   psd: power spectrum
%   freq: frequencies corresponding to the powers of the psd
%  Example usage: PLE = ple(ts, 500, [1 50], true)
% Authored by Yasir Çatal
% catalyasir@gmail.com
%% Do the calculation

% pwelch uses welch method that is useful when you have many sampling
% points (like EEG / MEG data) since it uses a sliding window approach when
% calculating PSD (take PSDs of chunks and average them, see the 
% documentation of pwelch for details). If your data has low number of time 
% points, you might want to switch to periodogram. A useful rule of thumb:
% your data should cover 2-3 full cycles of the lowest frequency you are
% interested in (for example, if your lowest frequency is 1 Hz, your chunks
% should be 2-3 seconds). 

%[psd, freq] = pwelch(ts, [], [], [], fs); 
[psd, freq] = periodogram(ts, [], [], fs);

% Next step is taking the frequencies you are interested in. The Nyquist
% for EEG data can go up to 500 Hz, but neuronally those frequencies are
% not thought to be very informative. 

foi = (freq > freqrange(1)) & (freq < freqrange(2));
freq2 = freq(foi);
psd2 = psd(foi);

% Finally, fit a straight line through the log-frequency log-power plot.
% Taking the log brings the alpha down from the exponent, turning the
% relation to a linear one (log(x^a) = a log(x)). Taking the negative is
% turns p = f^(alpha) to p = 1/f^(alpha)

coeffs = polyfit(log10(freq2), log10(psd2), 1);
PLE = -coeffs(1);

% Plot if you want
if isplot
    y = coeffs(2) + coeffs(1)*log10(freq2);
    loglog(freq, psd);
    hold on;
    loglog(freq2,10.^y,'r--');
    xlabel('Log Frequency')
    ylabel('Log Power')
    title(['PLE = ' num2str(PLE)])
    grid on
end
end