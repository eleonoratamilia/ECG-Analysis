%% Initialize TWA Analysis Pipeline
% This script sets up the necessary paths and processes ECG files for
% T-Wave Alternans analysis using spectral and MMA methods with HRV.

%% Add necessary libraries to MATLAB path
% Get the directory where THIS function is located
thisFile = mfilename('fullpath');
thisPath = fileparts(thisFile);

% Go one level up and add Helper Functions
parentPath = fileparts(thisPath);
helperPath = fullfile(parentPath, 'Helper Functions');
if exist(helperPath, 'dir')
    addpath(helperPath, '-begin');
    fprintf('✓ Helper Functions added to path\n');
else
    warning('Helper Functions folder not found at: %s', helperPath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add ECGdeli - HRVTool - Brainstorm to path%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Checking library dependencies...\n');

% Check if HRVTool is accessible
if ~isHRVToolAccessible()
    error('HRVTool not accessible. Please check HRVTool installation and path.');
else
    fprintf('✓ HRVTool is accessible\n');
end

% Check if ECGdeli is accessible
if ~isECGdeliAccessible()
    error('ECGdeli not accessible. Please check ECGdeli installation and path.');
else
    fprintf('✓ ECGdeli is accessible\n');
end

%% Display initialization status and user instructions
fprintf('TWA Analysis Pipeline initialized successfully.\n');
fprintf('Channel types must be named as ECG in Brainstorm.\n');
fprintf('• Supported labels: ECG, EKG (case-insensitive)\n');
fprintf('• Detection method: contains(lower(ch_type), {''ecg'', ''ekg''})\n\n');

%% Process ECG files for TWA analysis
fprintf('sFiles should be defined from Brainstorm: sFiles = { ... }\n\n');

% Check if sFiles exists and is not empty
if ~exist('sFiles', 'var') || isempty(sFiles)
    error('sFiles is not defined or empty. Please define sFiles from Brainstorm.');
end

% Process each file in the sFiles list
for iFile = 1:length(sFiles)
    fprintf('Processing file %d of %d: %s\n', iFile, length(sFiles), sFiles{iFile});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Step 1: Run TWA analysis pipeline   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [final_table, T_spectral_all, T_MMA_all, T_HRV_signal_all] = ...%%%%
        TWA_pipeline(sFiles{iFile}, '_twa_with_HRV'); %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Step 2: Extract study information for file naming
    [sStudy, iStudy, iDataFile, dataType, sItem] = bst_get('AnyFile', sFiles{iFile});
    filePath = sStudy.FileName;
    
    % Parse file path to extract subject/study identifier
    pathParts = strsplit(filePath, '/');
    if length(pathParts) >= 2
        rawIdentifier = pathParts{2};  % e.g., '@raw3a_EC_Sleep_Dec_2017'
        
        % Remove '@raw' prefix for clean file naming
        cleanIdentifier = regexprep(rawIdentifier, '@raw', '');
        
        %% Step 3: Save results with descriptive filename
        outputFilename = sprintf('%s_twa_with_HRV.mat', cleanIdentifier);
        save(outputFilename, 'final_table', 'T_spectral_all', 'T_MMA_all', 'T_HRV_signal_all');
        
        fprintf('Results saved to: %s\n\n', outputFilename);
    else
        warning('Unexpected file path structure: %s', filePath);
    end
end

%% Display completion message
fprintf('TWA analysis completed for all %d files.\n', length(sFiles));

%% Validation functions
function isAccessible = isHRVToolAccessible()
%ISHRVTOOLACCESSIBLE Check if HRVTool functions are available
    try
        % Try to access a core HRVTool function
        if exist('HRVTool', 'file') == 2
            isAccessible = true;
        else
            isAccessible = false;
        end
    catch
        isAccessible = false;
    end
end

function isAccessible = isECGdeliAccessible()
%ISECGDELIACCESSIBLE Check if ECGdeli functions are available
    try
        % Try to access a core ECGdeli function
        if exist('Annotate_ECG_Multi', 'file') == 2
            isAccessible = true;
        else
            isAccessible = false;
        end
    catch
        isAccessible = false;
    end
end