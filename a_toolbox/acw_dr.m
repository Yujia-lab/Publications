function [decayrate, acf, lags] = acw_dr(ts, fs, isplot, trim, trimtp)
%% Estimate the decay rate of ACF
% Input: 
%    ts : 1D vector 
%         Time series.
%    fs : double
%         Sampling rate (in Hz).
%    plot : Logical
%         If true, plots the fit to visually check
%    trim : Logical
%         If True, trim what comes after the ACF reaches 0.
%    trimtp : Integer
%         Number of time points to count after ACF reaches 0 before trimming. 
%         For example, if this is 5, the ACF will be trimmed after 5 lags from ACW-0
%         point. 
% Output: 
%    decayrate: Decay rate of ACF (in seconds)
%    acf:       Autocorrelation function
%    lags:      x-axis of ACF, for plotting purposes
arguments
    ts (:, 1)
    fs
    isplot = false;
    trim = true;
    trimtp = 5;
end

a = ver('econ');

if ~isempty(a) % Method using autocorr
    [acf, lags] = autocorr(ts, 'NumLags', length(ts)-1, 'NumSTD', 0);
else
    [acf, lags] = xcorr(ts, 'coeff');
    index = find(acf == 1); % Get rid of the left-side
    acf = acf(index:end); lags = lags(index:end);
end

lags = lags / fs; % Convert from samples to seconds
% Decoration: Get rid of the portion after ACW-0.
if trim
    [~, acw_0_i] = max(acf<=0);
    acf = acf(1:(acw_0_i + trimtp));
    lags = lags(1:(acw_0_i + trimtp));
end

% Initial guess (starting point) in nonlinear optimization. 
% Note that this is a bit arbitrary. What matters is that this has to be 
% same in all subjects
init_guess = 0.5;
expdecayfun = @(dr, x) exp(-dr*x);

decayrate = lsqcurvefit(expdecayfun, init_guess, lags, acf, 0, inf);

if isplot
    plot(lags, expdecayfun(decayrate, lags))
    hold on
    plot(lags, acf)
end

