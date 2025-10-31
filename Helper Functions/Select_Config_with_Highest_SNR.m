function Best_Config = Select_Config_with_Highest_SNR(Fs, signal)
%SELECT_CONFIG_WITH_HIGHEST_SNR Choose best ECG configuration based on SNR.
%   Best_Config = SELECT_CONFIG_WITH_HIGHEST_SNR(Fs, signal) evaluates multiple
%   signal configurations (single channels and bipolar combinations) and selects
%   the one with the highest Signal-to-Noise Ratio.
%
%   Inputs
%   ------
%   Fs     : positive scalar
%       Sampling frequency in Hz.
%   signal : numeric matrix
%       Multi-channel ECG signal (samples x channels).
%
%   Output
%   ------
%   Best_Config : numeric vector
%       Single-channel ECG signal with the highest SNR.

    %% Preprocess all channels
    Filtered_ECG_data = ECG_PreProcessing(Fs, signal);

    %% Handle single-channel case (early return)
    if size(signal, 2) == 1
        Best_Config = Filtered_ECG_data;
        return;
    end

    %% Evaluate SNR for each individual channel
    SNR_values = zeros(size(signal, 2), 1);
    for j = 1:size(Filtered_ECG_data, 2)
        SNR_values(j) = SNR_estim(Filtered_ECG_data(:, j), Fs);
    end
    
    % Find best single channel
    [best_SNR, best_channel] = max(SNR_values);
    Best_Config = Filtered_ECG_data(:, best_channel);

    %% Test bipolar configuration: raw subtraction then filter
    bipolar_raw_subtract = signal(:, 1) - signal(:, 2);
    bipolar_filtered = ECG_PreProcessing(Fs, bipolar_raw_subtract);
    bipolar_SNR = SNR_estim(bipolar_filtered, Fs);
    
    if bipolar_SNR > best_SNR
        best_SNR = bipolar_SNR;
        Best_Config = bipolar_filtered;
    end

    %% Test bipolar configuration: filter then subtract
    bipolar_filter_subtract = Filtered_ECG_data(:, 1) - Filtered_ECG_data(:, 2);
    bipolar_filter_subtract_SNR = SNR_estim(bipolar_filter_subtract, Fs);
    
    if bipolar_filter_subtract_SNR > best_SNR
        Best_Config = bipolar_filter_subtract;
    end
end