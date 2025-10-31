function [Values, index, start_indices, end_indices] = Spectral_twa(ecg, freq, Def_Params)
%SPECTRAL_TWA Calculate T-Wave Alternans using Spectral Method.
%   [Values, index, start_indices, end_indices] = SPECTRAL_TWA(ecg, freq, Def_Params)
%   computes T-wave alternans using spectral analysis on overlapping beat sequences.
%
%   Inputs
%   ------
%   ecg        : numeric vector
%       ECG signal in microvolts (uV)
%   freq       : positive scalar
%       Sampling frequency in Hz
%   Def_Params : struct
%       Parameter structure containing Spectral method configuration
%
%   Outputs
%   ------
%   Values        : numeric array (leads x methods x metrics x sequences)
%       TWA results for Lomb, Replace, and Difference methods
%   index         : numeric vector
%       QRS complex indices for each analyzed sequence
%   start_indices : numeric vector
%       Start sample indices for each TWA analysis window
%   end_indices   : numeric vector
%       End sample indices for each TWA analysis window

    %% Initialize global parameters
    clear global Param
    global Param beats lead
    
    lead = 1;
    Param.Metric = 'SM';  % Spectral Method
    Param.MethodForEctopy = Def_Params.Spectral.MethodForEctopy;
    Param.Alignment = 'st';
    beats = Def_Params.beats;

    %% Step 1: ECG annotation and preprocessing
    ann = Annotation_Table(ecg, freq);
    [s, tend, q, qs, stend] = Interpret_Annotations(ann, ecg);

    %% Step 2: Configure spectral method parameters
    Param.stAdjIntv = floor(Def_Params.Spectral.stAdjIntv_Multiplier * freq);
    Param.corrQRS = Def_Params.Spectral.corrQRS;      % QRS correlation threshold
    Param.corrT = Def_Params.Spectral.corrT;          % T-wave correlation threshold
    Param.RatioThreshold = Def_Params.Spectral.RatioThreshold;  % Significance threshold
    Param.incr = Def_Params.Spectral.incr;            % Beat increment between sequences
    Param.NInapr = Def_Params.Spectral.NInapr;        % Beat acceptance threshold

    %% Step 3: Validate beat count and prepare data
    % Require at least 90% of target beats for spectral analysis
    if length(s) < 0.9 * beats
        disp('Not enough QRS complexes for spectral analysis');
        Values = []; index = []; start_indices = []; end_indices = [];
        return;
    end
    
    % Replicate beats if needed to reach target count
    s = TWAUpdateArray(s, beats);
    q = TWAUpdateArray(q, beats);
    qs = TWAUpdateArray(qs, beats);
    stend = TWAUpdateArray(stend, beats);

    % Calculate ST segment length and apply filtering
    stlen = ApproximateSTLen(s, stend);
    ecg = FilterForTWA(ecg, freq);

    %% Step 4: Initialize spectral analysis variables
    clear global TWARes Align
    global TWARes Align
    Align.fid = [];
    StartQRSInd = 0;
    
    % Initialize result storage
    tracker = 1;
    Values = zeros(size(ecg, 2), 3, 2, tracker);  % leads x methods x metrics x sequences
    fid_points_saved = zeros(Def_Params.beats, tracker);
    QRSIND_tracker = [];
    start_indices = [];
    end_indices = [];

    %% Step 5: Spectral analysis on overlapping beat sequences
    while StartQRSInd + Param.incr + beats <= min(length(q), length(s))
        if StartQRSInd == 0
            StartQRSInd = 1;
        else
            StartQRSInd = StartQRSInd + Param.incr;
            disp(['TWASpectral: moving ' num2str(Param.incr) ' beats forward, starting from beat ' num2str(StartQRSInd)]);
        end

        % Align beats for current sequence
        disp('TWASpectral: aligning beats...');
        beat_range = StartQRSInd:StartQRSInd + beats - 1;
        Align = AlignBeats(ecg, beats, q(beat_range), s(beat_range), stlen);
        
        if isempty(Align.fid)
            disp('Alignment: failed');
            continue;
        else
            disp('Alignment: succeeded');
        end

        % Perform spectral TWA analysis
        disp('TWASpectral: looking for alternans...');
        [TWARes, Align] = TWASpectral_calculator(ecg, Align);
        
        if IsSuccessfull(TWARes)
            disp('TWASpectral: succeeded');
            
            % Store results for all three spectral methods
            method_names = {'lomb', 'replace', 'differences'};
            for i = 1:length(method_names)
                if ~isfield(TWARes, method_names{i}) || isempty(TWARes.(method_names{i}).significant)
                    continue;
                end
                Values(:, i, 1, tracker) = TWARes.(method_names{i}).VAlt;  % Alternans magnitude
                Values(:, i, 2, tracker) = TWARes.(method_names{i}).Ratio; % Significance ratio
            end
            
            % Store sequence information
            fid_points_saved(:, tracker) = Align.fid(:);
            tracker = tracker + 1;
            QRSIND_tracker = [QRSIND_tracker; StartQRSInd];
            start_indices = [start_indices; Align.fid(1)];
            end_indices = [end_indices; Align.fid(end)];
        end
    end

    %% Step 6: Format final outputs
    index = s(QRSIND_tracker);  % QRS complex indices for each sequence
end

function res = IsSuccessfull(TWARes)
%ISSUCCESSFULL Check if spectral TWA analysis produced valid results.
%   Returns true if any of the three methods (Lomb, Replace, Differences)
%   successfully detected T-wave alternans.

    res = (isfield(TWARes, 'lomb') && TWARes.lomb.successfull) || ...
          (isfield(TWARes, 'replace') && TWARes.replace.successfull) || ...
          (isfield(TWARes, 'differences') && TWARes.differences.successfull);
end