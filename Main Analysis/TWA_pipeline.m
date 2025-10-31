function [final_table, T_spectral_all, T_MMA_all, T_HRV_signal_all] = TWA_pipeline(File_name, varargin)
%TWAPIPELINE Complete T-Wave Alternans analysis pipeline for Brainstorm data.
%   [final_table, T_spectral_all, T_MMA_all, T_HRV_signal_all] = TWAPIPELINE(File_name)
%   performs comprehensive TWA analysis including spectral method, MMA method,
%   and HRV analysis on Brainstorm ECG data.
%
%   Inputs
%   ------
%   File_name : string
%       Brainstorm record name to analyze
%   varargin  : optional extension suffix for output files
%
%   Outputs
%   ------
%   final_table       : Summary table with aggregated TWA and HRV metrics
%   T_spectral_all    : Detailed spectral TWA results table
%   T_MMA_all         : Detailed MMA TWA results table  
%   T_HRV_signal_all  : Detailed HRV analysis results table

    %% Initialize file extension
    if length(varargin) >= 1
        extension = varargin{1};
    else
        extension = '_TWA_res';
    end
    
    %% Get record details from Brainstorm database
    [sStudy, iStudy, iFile, DataType, sItem] = bst_get('AnyFile', File_name);
    [sSubject, iSubject] = bst_get('Subject', sStudy.BrainStormSubject, 1);

    %% Import and split data for memory-efficient processing
    sFiles_imported = bst_process('CallProcess', 'process_import_data_time', File_name, [], ...
        'subjectname',   sSubject.Name, ...
        'condition',     'TWA3', ...
        'timewindow',    [], ...
        'split',         3600, ... % Split into 1-hour segments for memory management
        'ignoreshort',   0);

    [~, iStudy, ~, ~, ~] = bst_get('AnyFile', sFiles_imported(1).FileName);

    %% Define TWA analysis channels
    % Get channel file
    channelFile = bst_get('ChannelFileForStudy', File_name);
    channelInfo = in_bst_channel(channelFile);
    nExisting = length(channelInfo.Channel);

    % Define new channels
    twaChannels = {
        'Spec_Max',        'TWA',          'Spectral TWA maximum amplitude'
        'MMA',             'TWA',          'MMA TWA amplitude' 
        'E_C_G_filtered',  'E_C_G_filtered', 'Filtered ECG signal'
        'TWA_Spectral_Ratio', 'TWA',       'Spectral TWA significance ratio'
        'HR',              'HR',           'Heart rate'
        'HRV_SDNN',        'HRV',          'HRV SDNN metric'
        'HRV_RMSSD',       'HRV',          'HRV RMSSD metric'
        'HRV_rrHRV',       'HRV',          'HRV rrHRV metric'
    };
    
    % Create new channels
    for iChan = 1:size(twaChannels, 1)
        channelInfo.Channel(nExisting + iChan).Name   = twaChannels{iChan, 1};
        channelInfo.Channel(nExisting + iChan).Type   = twaChannels{iChan, 2};
        channelInfo.Channel(nExisting + iChan).Loc    = [0; 0; 0];
        channelInfo.Channel(nExisting + iChan).Weight = 1;
    end

    % Save channel file - FIXED: changed 'ch' to 'channelInfo'
    db_set_channel(iStudy, channelInfo, 2, 0);

    %% Detect ECG channels
    channelTypes = {channelInfo.Channel.Type};
    ecgIndices = find(contains(lower(channelTypes), {'ecg', 'ekg'}));
    
    %% Case with no ECG channels: Return Empty
    if isempty(ecgIndices)
        fprintf('No ECG channels detected. Returning empty results.\n');
        [final_table, T_spectral_all, T_MMA_all, T_HRV_signal_all] = createEmptyResults();
        initializeEmptySignalFiles(sFiles_imported, nExisting);
        return;
    end

    %% Initialize result storage
    T_MMA_all = [];
    T_spectral_all = [];
    T_HRV_signal_all = [];
    RR_all = [];
    History_save = [];

    %% Process each imported file segment
    for i = 1:length(sFiles_imported)
        file_name_imp = sFiles_imported(i).FileName;
        file_path_imp = file_fullpath(file_name_imp);
        
        % Load segment data
        segmentData = load(file_path_imp);
        Time = segmentData.Time;
        F = segmentData.F;

        % Calculate sampling frequency
        Fs = calculateSamplingFrequency(Time);

        % Select the ECG configuration with highest SNR
        ecgData = F(ecgIndices, :)' * 1e6;  % Convert to microvolts
        bestConfig = Select_Config_with_Highest_SNR(Fs, ecgData);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Perform TWA analysis on current segment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [T_spectral, T_MMA, T_HRV_signal, RR] = TWA_process(bestConfig, Time, Fs, i);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Aggregate results
        T_MMA_all = [T_MMA_all; T_MMA];
        T_spectral_all = [T_spectral_all; T_spectral];
        T_HRV_signal_all = [T_HRV_signal_all; T_HRV_signal];
        RR_all = [RR_all; RR];

        %% Write TWA results to signal data
        F = writeTWAResultsToSignal(F, nExisting, T_spectral, T_MMA, bestConfig, T_HRV_signal);

        % Update channel flags and history
        segmentData.ChannelFlag = [segmentData.ChannelFlag; ones(8, 1)];
        segmentData.F = F; % Update the signal data

        if i == 1
            History = segmentData.History;
            a = size(History, 1);
            History{a + 1, 1} = History{a, 1};
            History{a + 1, 2} = 'Do_TWA_Analysis_Ege';
            History{a + 1, 3} = 'Updated Manually';
            History_save = History;
        end
        
        % Save updated segment data
        save(file_path_imp, '-struct', 'segmentData', ...
             'ChannelFlag', 'ColormapType', 'Comment', 'DataType', 'Device', ...
             'DisplayUnits', 'Events', 'F', 'History', 'Leff', 'nAvg', 'Std', 'Time');
    end

    %% Reconstruct and finalize results
    [~, iStudy, ~, ~, sItem] = bst_get('AnyFile', file_name_imp);
    db_reload_studies(iStudy);
    
    sFiles = [];
    for i = 1:length(sFiles_imported)
        file_name_imp = sFiles_imported(i).FileName;
        % Transform Imported file into Raw File
        [sStudy, iStudy, iFile, DataType, sItem] = bst_get('AnyFile', file_name_imp);
        [sSubject, iSubject] = bst_get('Subject', sStudy.BrainStormSubject, 1);
    
        new_file_name = import_raw(file_fullpath(file_name_imp), 'BST-DATA', iSubject);
        sFiles{1, i} = new_file_name{1}; % The imported raw ones
    end
    
    sFiles_final = bst_process('CallProcess', 'process_concat', sFiles, []);

    %% Update the History
    file_path = file_fullpath(sFiles_final.FileName);
    data = in_bst_data(sFiles_final.FileName);
    History = History_save;
    save(file_path, 'History', '-append');
    %% Upload the Events file
    if ispc
        % Windows OS
        fileSeparator = '\';
    elseif isunix
        % Linux or macOS
        fileSeparator = '/';
    else
        % Fallback 
        fileSeparator = '_';
    end
    
    for i = 1:length(sFiles_imported)
        save_name = strcat(pwd, fileSeparator, 'events_', num2str(i), '.mat');

        bst_process('CallProcess', 'process_evt_import', sFiles_final.FileName, [], ...
            'evtfile', {save_name, 'BST'}, ...
            'evtname', '', ...
            'delete',  0);

        delete(save_name);
    end
    %% Change the name of the final study
    % Get the new study name
    division = regexp(File_name, '[\\/]');
    study_name = File_name(division(1) + 1:division(2) - 1);
    if startsWith(study_name, '@raw')
        study_name = strcat(study_name, extension);
    else
        study_name = strcat('@raw', study_name, extension);
    end
    
    % Reload the study
    [~, iStudy, ~, ~, sItem] = bst_get('AnyFile', sFiles_final.FileName);
    db_reload_studies(iStudy);
    
    % Change the name
    old_path = sFiles_final.FileName(1:max(regexp(sFiles_final.FileName, '[\\/]')) - 1);
    new_path_all = old_path(1:min(regexp(old_path, '[\\/]')));
    new_path_all = strcat(new_path_all, study_name);
    db_rename_condition(old_path, new_path_all, 1, 1);
    %% Delete partitioned subfiles
    sFiles = bst_process('CallProcess', 'process_delete', sFiles, [], 'target', 2); 
    sFiles = bst_process('CallProcess', 'process_delete', {sFiles_imported(end).FileName}, [], 'target', 2);

    %% Generate final summary table
    final_table = Generate_TWA_HRV_SummaryTable(length(ecgIndices), T_spectral_all, T_MMA_all, RR_all);
end

%% Helper functions
function [final_table, T_spectral_all, T_MMA_all, T_HRV_signal_all] = createEmptyResults()
%CREATEEMPTYRESULTS Create empty result structures for no-ECG case
    final_table = table('Size', [0, 9], 'VariableTypes', repmat({'double'}, 1, 9), ...
        'VariableNames', {'Ecg_count', 'Spec_Max', 'MMA_Max', 'Spec_Median', 'MMA_Median', ...
                          'HR', 'HRV_SDNN', 'HRV_RMSSD', 'HRV_rrHRV'});
    
    final_table.Ecg_count(1) = 0;
    final_table.Spec_Max(1) = NaN;
    final_table.MMA_Max(1) = NaN;
    final_table.Spec_Median(1) = NaN;
    final_table.MMA_Median(1) = NaN;
    final_table.HR = NaN;
    final_table.HRV_SDNN = NaN;
    final_table.HRV_RMSSD = NaN;
    final_table.HRV_rrHRV = NaN;
    
    T_spectral_all = table();
    T_MMA_all = table();
    T_HRV_signal_all = table();
end

function initializeEmptySignalFiles(sFiles_imported, nExisting)
%INITIALIZEEMPTYSIGNALFILES Initialize signal files with placeholder values
    for i = 1:length(sFiles_imported)
        file_name_imp = sFiles_imported(i).FileName;
        file_path_imp = file_fullpath(file_name_imp);
        segmentData = load(file_path_imp);
        
        segmentData.F(nExisting + 1:nExisting + 8, :) = -1 * ones(8, size(segmentData.F, 2));
        segmentData.ChannelFlag = [segmentData.ChannelFlag; ones(8, 1)];
        
        save(file_path_imp, '-struct', 'segmentData', ...
             'ChannelFlag', 'ColormapType', 'Comment', 'DataType', 'Device', ...
             'DisplayUnits', 'Events', 'F', 'History', 'Leff', 'nAvg', 'Std', 'Time');
    end
end

function Fs = calculateSamplingFrequency(Time)
%CALCULATESAMPLINGFREQUENCY Calculate sampling frequency from time vector
    fsReal = (length(Time) - 1) / (Time(end) - Time(1));
    fsInt = round(fsReal);
    
    if abs(fsReal - fsInt) < max(1e-6 * fsInt, 1e-9)
        Fs = fsInt;   % Return clean integer
    else
        Fs = fsReal;  % Handle non-integer case
    end
end

function F = writeTWAResultsToSignal(F, nExisting, T_spectral, T_MMA, bestConfig, T_HRV_signal)
%WRITETWARESULTSTOSIGNAL Write TWA analysis results to signal matrix
    % Initialize TWA channels
    for iChan = 1:8
        F(nExisting + iChan, :) = zeros(1, size(F, 2));
    end
    
    % Write spectral TWA results
    if ~isempty(T_spectral)
        for j = 1:height(T_spectral)
            if isnan(T_spectral.Value(j, 1))
                F(nExisting + 1, T_spectral.Indices(j, 1):T_spectral.Indices(j, 2)) = -1;
            else
                F(nExisting + 1, T_spectral.Indices(j, 1):T_spectral.Indices(j, 2)) = T_spectral.Value(j, 1);
            end
            F(nExisting + 4, T_spectral.Indices(j, 1):T_spectral.Indices(j, 2)) = T_spectral.Ratio(j, 1);
        end
    end
    
    % Write MMA TWA results
    if ~isempty(T_MMA)
        for j = 1:height(T_MMA)
            F(nExisting + 2, T_MMA.Indices(j, 1):T_MMA.Indices(j, 2)) = T_MMA.Value(j, 1);
        end
    end
    
    % Write filtered ECG
    F(nExisting + 3, :) = bestConfig;
    
    % Write HRV results with linear interpolation
    if ~isempty(T_HRV_signal) && length(T_HRV_signal.Indices) >= 2
        for j = 2:length(T_HRV_signal.Indices)
            idxRange = T_HRV_signal.Indices(j-1):T_HRV_signal.Indices(j);
            nPoints = length(idxRange);
            
            F(nExisting + 5, idxRange) = linspace(T_HRV_signal.HR(j-1), T_HRV_signal.HR(j), nPoints);
            F(nExisting + 6, idxRange) = linspace(T_HRV_signal.SDNN(j-1), T_HRV_signal.SDNN(j), nPoints);
            F(nExisting + 7, idxRange) = linspace(T_HRV_signal.RMSSD(j-1), T_HRV_signal.RMSSD(j), nPoints);
            F(nExisting + 8, idxRange) = linspace(T_HRV_signal.rrHRV(j-1), T_HRV_signal.rrHRV(j), nPoints);
        end
    end
end