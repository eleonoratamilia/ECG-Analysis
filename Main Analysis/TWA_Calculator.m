function [final_table,T_spectral,T_MMA,T_bad] = TWA_Calculator(ecg,time,Fs)
%TWACALCULATOR Perform comprehensive T-Wave Alternans analysis on ECG data.
%   [final_table, T_spectral, T_MMA, T_bad] = TWACALCULATOR(ecg, time, Fs)
%   computes TWA using spectral and MMA methods with signal quality assessment.
%   also includes HRV information in final_table
%
%   Inputs
%   ------
%   ecg  : numeric matrix (time_points x channels)
%       ECG signal data in microvolts (uV)
%   time : numeric vector
%       Time vector corresponding to ECG samples
%   Fs   : positive scalar
%       Sampling frequency in Hz
%
%   Outputs
%   ------
%   final_table : table
%       Summary table with TWA and HRV metrics
%   T_spectral  : table
%       Detailed spectral TWA results
%   T_MMA       : table
%       Detailed MMA TWA results
%   T_bad       : table
%       Bad quality segments information

    %% Parameter Definition for TWA Analysis
    % Thresholds and configuration
    Def_Params.Spectral.corrQRS = 0.95;              % QRS correlation threshold
    Def_Params.Spectral.corrT = 0.8;                 % T-wave correlation threshold
    Def_Params.Spectral.RatioThreshold = 0;          % Significance threshold (applied later)
    Def_Params.Spectral.incr = 20;                   % Beat increment for spectral analysis
    Def_Params.Spectral.NInapr = 0.1;                % Beat acceptance threshold
    Def_Params.Spectral.MethodForEctopy = 'lomb replace differences'; % Spectral methods
    Def_Params.Spectral.stAdjIntv_Multiplier = 0.03; % Beat alignment tolerance
    
    Def_Params.beats = 128;                          % Beats for PSD calculation
    
    % MMA parameters
    Def_Params.MMA.stAdjIntv_Multiplier = 0.04;      % Beat alignment tolerance
    Def_Params.MMA.Interval = 15;                    % MMA result interval (seconds)
    Def_Params.MMA.corrQRS = 0.95;                   % QRS correlation threshold
    Def_Params.MMA.corrT = 0.8;                      % T-wave correlation threshold
    Def_Params.MMA.NInapr = 0.1;                     % Beat acceptance threshold
    
    % Signal segmentation parameters
    Def_Params.Spectral.minGoodLength = 120;         % Minimum good segment (seconds)
    Def_Params.Spectral.splitLength = 240;           % Target segment length (seconds)
    Def_Params.Spectral.maxGoodLength = 480;         % Maximum segment length (seconds)
    
    %% Initialize result storage
    T_spectral = table();
    T_MMA = table();
    T_HRV = table();
    T_bad = table();
    final_table =  table();
    time = time(:) ; % Ensure time is column vector

    %% Channel Count Validation
    if size(ecg, 2) > 20
        fprintf('Warning: More than 20 channels detected.\n');
        fprintf('Please verify this is intended.\n\n');
    end
    
    %% Select the ECG configuration with highest SNR
    ecg_filt = Select_Config_with_Highest_SNR(Fs, ecg);
    %% Signal Quality Check    
    [flags_res,good_quality] = sig_quality(ecg_filt,Fs,10,0.2);
    % Identify bad quality segments
    bad_quality = not(good_quality);
    changes = diff([0; bad_quality; 0]);  % Detect changes from 0 to 1 and 1 to 0
    startIndices = find(changes == 1);
    endIndices = find(changes == -1) - 1;
    T_bad.bad_times_start = time(startIndices);
    T_bad.bad_times_end = time(endIndices);
    
    [segmentationTable] = Split_signal(good_quality,Fs,Def_Params);

    %% RR for HRV Calculation
    RR=[];
    for i=1:size(segmentationTable,1)   
    ecg_samp = ecg_filt(segmentationTable{i,2} : segmentationTable{i,3}); % The good signal part
    [~,FPT_Cell_1]=Annotate_ECG_Multi(ecg_samp,Fs); % Delineate the ECG
    FPT_1 = FPT_Cell_1{1,1}; % Get the Fiducial Point Table
    R_peaks = time(FPT_1(:,6))';
    RR = [RR;diff(R_peaks(:))];
    end

    %% Spectral TWA Calculation
    warning('off', 'MATLAB:table:RowsAddedExistingVars');
    for i = 1:size(segmentationTable, 1)   
        segment_start = segmentationTable{i, 2};
        segment_end = segmentationTable{i, 3};
        ecg_segment = ecg_filt(segment_start:segment_end);
        shift_required = segment_start - 1;

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        [Values, index, start_indices, end_indices] = Spectral_twa(ecg_segment, Fs, Def_Params);
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        % Update table with spectral TWA results if alignment successful
        if ~isempty(index)
            start_indices = start_indices + shift_required;
            end_indices = end_indices + shift_required;
            
            tracker = height(T_spectral);
            for j = 1:length(index)
                [Y, I] = max(Values(1, :, 1, j));
                
                T_spectral.Indices(tracker + j, 1:2) = [start_indices(j), end_indices(j)];
                T_spectral.Time(tracker + j, 1:2) = [time(start_indices(j)), time(end_indices(j))];
                
                if Y == -1
                    T_spectral.Value(tracker + j, 1) = NaN;
                else
                    T_spectral.Value(tracker + j, 1) = Y;
                end
                
                T_spectral.Ratio(tracker + j, 1) = Values(1, I, 2, j);
            end
        end
    end

    %% MMA TWA Calculation 
    for i = 1:size(segmentationTable, 1)   
        segment_start = segmentationTable{i, 2};
        segment_end = segmentationTable{i, 3};
        ecg_segment = ecg_filt(segment_start:segment_end);
        shift_required = segment_start - 1;
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        [TWARes, index, end_indices] = MMA_twa(ecg_segment, Fs, Def_Params);
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        % Update table with MMA TWA results if alignment successful
        tracker = height(T_MMA);
        if isfield(TWARes, 'res')
            index = index + shift_required;
            end_indices = end_indices + shift_required;
            
            T_MMA.Indices(tracker + 1:tracker + size(TWARes.VAltTrend, 1), 1:2) = [index(:), end_indices(:)];
            T_MMA.Time(tracker + 1:tracker + size(TWARes.VAltTrend, 1), 1:2) = [time(index(:)), time(end_indices(:))];
            T_MMA.Value(tracker + 1:tracker + size(TWARes.VAltTrend, 1), 1) = max(TWARes.VAltTrend, [], 2);          
        end       
    end

    %% Generate final summary table
    final_table = Generate_TWA_HRV_SummaryTable(size(ecg,2), T_spectral, T_MMA, RR);
end