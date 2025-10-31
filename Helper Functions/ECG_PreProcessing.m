function [ecg_filtered_isoline] = ECG_PreProcessing(Fs, signal)
%ECG_PREPROCESSING Apply complete ECG preprocessing pipeline.
%   [ecg_filtered_isoline] = ECG_PREPROCESSING(Fs, signal) processes raw ECG
%   through baseline removal, notch filtering, bandpass filtering, and isoline
%   correction.
%
%   Requires ECGDeli functions in path.
%
%   Inputs
%   ------
%   Fs : positive scalar
%       Sampling frequency in Hz.
%   signal : numeric vector
%       Raw ECG signal to process.
%
%   Output
%   ------
%   ecg_filtered_isoline : numeric vector
%       Fully preprocessed ECG signal ready for analysis.

    %% Step 1: Baseline wander removal using median filtering
    [ecg_filt_baseline, ~] = ECG_Baseline_Removal(signal, Fs, 1, 0.5);
    
    %% Step 2: Powerline interference removal (60 Hz and harmonics)
    ecg_filtered_frq = Notch_Filter(ecg_filt_baseline, Fs, 60, 2);  % 60 Hz notch
    ecg_filtered_frq = Notch_Filter(ecg_filtered_frq, Fs, 120, 2);  % 120 Hz notch
    
    %% Step 3: Bandpass filtering for ECG frequency range
    [ecg_filtered_frq] = ECG_High_Low_Filter(ecg_filtered_frq, Fs, 1, 40);
    
    %% Step 4: Isoline correction (subtract median of entire channel)
    [ecg_filtered_isoline, ~, ~, ~] = Isoline_Correction(ecg_filtered_frq);
end