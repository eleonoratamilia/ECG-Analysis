function [flags_res, good_quality] = sig_quality(ecg, freq, Step_time, flatline_time)
%SIG_QUALITY Assess ECG signal quality across sliding windows.
%   [flags_res, good_quality] = SIG_QUALITY(ecg, freq, Step_time, flatline_time)
%   evaluates signal quality based on flatline detection, SNR, and heart rate
%   criteria across sequential time windows.
%
%   Inputs
%   ------
%   ecg          : numeric vector
%       ECG signal to analyze
%   freq         : positive scalar
%       Sampling frequency in Hz
%   Step_time    : positive scalar
%       Window duration for analysis in seconds
%   flatline_time : positive scalar
%       Sub-window duration for flatline detection in seconds
%
%   Outputs
%   ------
%   flags_res    : logical vector
%       Quality flags for each window (true = good quality)
%   good_quality : numeric vector
%       Binary mask indicating good quality samples (1 = good, 0 = poor)

    %% Initialize variables
    good_quality = zeros(length(ecg), 1);
    duration = floor(Step_time * freq);
    
    %% Detect R-peaks for heart rate calculation
    fpt = QRS_Detection(ecg, freq, 'peaksQRS', 'm');
    R_wave = fpt(:, 6);

    %% Analyze signal quality in sequential windows
    n_step = floor(length(ecg) / duration);
    flags = false(3, n_step);  % [flatline; SNR; heart_rate]
    
    for i = 1:n_step
        start_sample = (i - 1) * duration + 1;
        end_sample = i * duration;
        ecg_window = ecg(start_sample:end_sample);
        
        % Criterion 1: Flatline detection
        bb = buffer_with_truncate(ecg_window, round(freq * flatline_time));
        flags(1, i) = ~any(arrayfun(@(x) check_flatline(bb(:, x)), 1:size(bb, 2)));
        
        % Criterion 2: Signal-to-Noise Ratio
        % Threshold 1.122 corresponds to 1 dB power ratio
        flags(2, i) = SNR_estimate(ecg_window, freq) > 1.122;
        
        % Criterion 3: Heart rate validity (24-300 BPM)
        R_in_window = R_wave(R_wave >= start_sample & R_wave <= end_sample);
        
        if length(R_in_window) < 2
            flags(3, i) = false;
            continue;
        end
        
        RR_intervals = diff(R_in_window);
        HR = 60 * freq / mean(RR_intervals);
        flags(3, i) = (HR >= 24 && HR <= 300);
        
        % Mark good quality segments
        if all(flags(:, i))
            good_quality(start_sample:end_sample) = 1;
        end
    end
    
    flags_res = all(flags, 1);
end