function SNR = SNR_estimate(sig, Fs)
%SNR_estimate Estimate Signal-to-Noise Ratio using power spectral density.
%   SNR = SNR_ESTIM(sig, Fs) calculates SNR by comparing signal power in the
%   ECG frequency band (2-40 Hz) against noise power in low and high frequencies.
%
%   Inputs
%   ------
%   sig : numeric vector
%       Input ECG signal
%   Fs  : positive scalar
%       Sampling frequency in Hz
%
%   Output
%   ------
%   SNR : positive scalar
%       Signal-to-Noise Ratio (signal_power / noise_power)

    %% Calculate power spectral density using periodogram
    [pxx, f] = periodogram(sig, [], [], Fs);
    
    %% Define frequency bands for signal and noise
    % Signal band: 2-40 Hz (typical ECG frequencies)
    % Noise bands: 0-2 Hz (baseline wander) and 40-250 Hz (high-frequency noise)
    
    idx_0Hz   = find(f >= 0,   1, 'first');   % Start of spectrum
    idx_2Hz   = find(f >= 2,   1, 'first');   % Start of signal band
    idx_40Hz  = find(f >= 40,  1, 'first');   % End of signal band
    idx_250Hz = find(f >= 250, 1, 'first');   % End of noise band
    
    %% Calculate signal and noise power
    signal_power = sum(pxx(idx_2Hz:idx_40Hz));
    noise_power  = sum(pxx(idx_0Hz:idx_2Hz)) + sum(pxx(idx_40Hz:idx_250Hz));
    
    %% Compute Signal-to-Noise Ratio
    SNR = signal_power / noise_power;
end