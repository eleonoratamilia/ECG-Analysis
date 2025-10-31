function [T_spectral, T_MMA, T_HRV_signal, RR] = TWA_process(ecg_filt, time, Fs, i)
%TWAPROCESS Perform T-Wave Alternans analysis on ECG signal.
%   [T_spectral, T_MMA, T_HRV_signal, RR] = TWAPROCESS(ecg_filt, time, Fs, i)
%   computes TWA using spectral and MMA methods with HRV analysis.
%
%   Inputs
%   ------
%   ecg_filt : filtered ECG signal
%   time     : time vector
%   Fs       : sampling frequency
%   i        : file identifier
%
%   Outputs
%   ------
%   T_spectral    : Spectral TWA results
%   T_MMA         : MMA TWA results
%   T_HRV_signal  : HRV analysis results
%   RR            : RR intervals

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

    %% Signal Quality Assessment
    [flags_res, good_quality] = sig_quality(ecg_filt, Fs, 10, 0.2);
    [segmentationTable] = Split_signal(good_quality, Fs, Def_Params);

    %% Create All ECG Analysis Events
    warning('off', 'MATLAB:table:RowsAddedExistingVars');
    events = createECGEvents(ecg_filt, time, Fs, i, good_quality);

    %% Calculate HRV signal
    [~,FPT_Cell]=Annotate_ECG_Multi(ecg_filt,Fs); % Delineate the ECG
    FPT = FPT_Cell{1,1}; % Get the Fiducial Point Table
    R_peaks = time(FPT(:,6))';
    R_indices = FPT(:,6);
    
    RR_whole_sig = [NaN;diff(R_peaks(:))];
    RR_whole_sig = HRV.RRfilter(RR_whole_sig,0.15);

    HR = HRV.HR(RR_whole_sig,60);    
    SDNN = HRV.SDNN(RR_whole_sig,60);   
    RMSSD = HRV.RMSSD(RR_whole_sig,60);  
    rrHRV = HRV.rrHRV(RR_whole_sig,60);
    
    HR(isnan(HR)) = 0;
    SDNN(isnan(SDNN)) = 0;
    RMSSD(isnan(RMSSD)) = 0;
    rrHRV(isnan(rrHRV)) = 0;

    T_HRV_signal = table(R_indices(:),R_peaks(:),HR,SDNN,RMSSD,rrHRV,'VariableNames',{'Indices','Time','HR','SDNN','RMSSD','rrHRV'});

    %% Calculate RR Intervals from Good Quality Segments Only
    RR = [];
    for i=1:size(segmentationTable,1)   
        segment_start = segmentationTable{i, 2};
        segment_end = segmentationTable{i, 3};
        ecg_segment = ecg_filt(segment_start:segment_end);
        
        [~, FPT_Cell_1] = Annotate_ECG_Multi(ecg_segment, Fs);
        FPT_1 = FPT_Cell_1{1, 1};
        R_peaks = time(FPT_1(:, 6))';
        RR = [RR; diff(R_peaks(:))];
    end

    %% Making Sure time is a column vector
    time = time(:);

   %% Spectral TWA Calculation
    T_spectral = table(); % Initialize results table
    
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
    T_MMA = table(); % Initialize results table
    
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
            
            T_MMA.Indices(tracker + 1:tracker + size(TWARes.VAltTrend, 1), 1:2) = ...
                [index(:), end_indices(:)];
            T_MMA.Time(tracker + 1:tracker + size(TWARes.VAltTrend, 1), 1:2) = ...
                [time(index(:)), time(end_indices(:))];
            T_MMA.Value(tracker + 1:tracker + size(TWARes.VAltTrend, 1), 1) = ...
                max(TWARes.VAltTrend, [], 2);          
        end       
    end
end


function events = createECGEvents(ecg_filt, time, Fs, fileID, good_quality)
%CREATEECGEVENTS Create comprehensive ECG analysis events structure
%   Creates events for signal quality, R-peaks, and waveform delineation

    %% Initialize events structure
    events = [];
    
    %% Create file path
    if ispc
        fileSeparator = '\';
    elseif isunix
        fileSeparator = '/';
    else
        fileSeparator = '_';
    end
    save_name = strcat(pwd, fileSeparator, 'events_', num2str(fileID), '.mat');

    %% Find Bad Signal Events
    bad_quality = not(good_quality);
    changes = diff([0; bad_quality; 0]);             % Detect quality transitions
    startIndices = find(changes == 1);
    endIndices = find(changes == -1) - 1;

    %% 1. Signal Quality Events
    tim(1, 1:length(startIndices)) = time(startIndices);
    tim(2, 1:length(startIndices)) = time(endIndices);

    events(1).label = 'ECG_sig_quality';
    events(1).color = [0.0200, 0.0200, 1];
    events(1).epochs = ones(1, length(startIndices));
    events(1).times = tim;
    events(1).reactTimes = [];
    events(1).select = 1;
    events(1).channels = []; 
    events(1).notes = [];

    %% 2. ECG Delineation and Waveform Events
    [~, FPT_Cell] = Annotate_ECG_Multi(ecg_filt, Fs);
    FPT = FPT_Cell{1, 1};

    % Clean Fiducial Point Table
    n = length(ecg_filt);
    if isempty(FPT(end, 12)) || FPT(end, 12) > n
        FPT(end, :) = []; % Remove incomplete last beat
    end
    if isempty(FPT(1, 1)) || FPT(1, 1) < 2
        FPT(1, :) = []; % Remove first beat if P-wave missing
    end
    rowsWithZero = (any(FPT(:, 1:12) == 0, 2));
    FPT(rowsWithZero, :) = []; % Remove failed detections

    % Define waveform configurations [label, color, FPT_start_col, FPT_end_col]
    waveforms = {
        'R_peak',  [0.56, 0.01, 0.91], 6, 6    % Single point
        'P_wave',  [1, 0, 1],           1, 3    % Start to peak
        'QRS_wave',[0, 1, 0],           4, 8    % Onset to offset  
        'T_wave',  [0, 1, 1],           10, 12  % Peak to end
    };

    % Create waveform events
    for iWave = 1:size(waveforms, 1)
        eventIndex = iWave + 1; % +1 because event(1) is signal quality
        
        events(eventIndex).label = waveforms{iWave, 1};
        events(eventIndex).color = waveforms{iWave, 2};
        events(eventIndex).epochs = ones(1, size(FPT, 1));
        
        if waveforms{iWave, 3} == waveforms{iWave, 4}
            % Single point event (R-peak)
            events(eventIndex).times = time(FPT(:, waveforms{iWave, 3}))';
        else
            % Range event (P-wave, QRS, T-wave)
            events(eventIndex).times = [time(FPT(:, waveforms{iWave, 3}))'; ...
                                       time(FPT(:, waveforms{iWave, 4}))'];
        end
        
        events(eventIndex).reacTimes = [];
        events(eventIndex).select = 0;
        events(eventIndex).channels = [];
        events(eventIndex).notes = [];
    end

    %% Save events to file
    save(save_name, 'events');
end