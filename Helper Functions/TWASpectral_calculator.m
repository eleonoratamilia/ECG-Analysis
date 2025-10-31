function [TWARes, Align] = TWASpectral_calculator(ecg, Align)
%TWASPECTRAL_CALCULATOR Compute T-Wave Alternans using Spectral Method.
%   [TWARes, Align] = TWASPECTRAL_CALCULATOR(ecg, Align) calculates T-wave
%   alternans using spectral analysis on aligned ECG beats with multiple
%   methods for handling ectopic beats and noise.
%
%   Inputs
%   ------
%   ecg   : numeric matrix (samples x leads)
%       Multi-lead ECG signal data
%   Align : struct
%       Beat alignment structure from AlignBeats containing:
%         - fid:      Fiducial points (QRS positions)
%         - f2s:      Offset to ST segment start
%         - amp:      Beat amplitudes for baseline correction
%         - st:       ST segment length in samples
%         - valid:    Beat validity flags
%         - validleads: Lead validity flags
%
%   Outputs
%   ------
%   TWARes : struct
%       Results structure containing TWA measurements with fields for each method:
%         - lomb:      Lomb periodogram method for unevenly sampled data
%         - replace:   Ectopic beat replacement method
%         - differences: Beat-to-beat difference method (natural detrending)

% TWARes: structure containing the fields:
%       alt_series     TWA sequence array
%       avg_even (num_of_timepoints x num_of_leads): time average of even beats
%       avg_odd  (num_of_timepoints x num_of_leads): time average of odd beats
%       psd      (num_of_beats/2+1 x num_of_leads x num_of_timepoints): psd of timepoints within Q-Tend
%       lomb     (num_of_beats/2+1 x num_of_leads x num_of_timepoints): lomb psd of timepoints within Q-Tend
%       avg_psd  (num_of_beats/2+1 x num_of_leads): average psd across all timepoints in Q-Tend
%       significant    (1 x num_of_leads): 1 if the alternans value is statistically significant against the Param.RatioThreshold in the lead and 0 otherwise
%       VAlt     (1 x num_of_leads):
%       Ratio    (1 x num_of_leads):


    global Param
    TWARes = [];

    %% Initialize result structures for requested PSD methods
    methodNames = {'lomb', 'replace', 'differences'};
    for i = 1:length(methodNames)
        if contains(Param.MethodForEctopy, methodNames{i}, 'IgnoreCase', true)
            TWARes.(methodNames{i}).successfull = false;
        end
    end
   
    %% Process each ECG lead
    for lead = 1:size(ecg, 2)
        % Skip leads marked as invalid during alignment
        if ~Align.validleads(lead)
            continue;
        end
    
     %% Method 1: Lomb Periodogram (handles uneven beat intervals)
        if contains(Param.MethodForEctopy, 'lomb','IgnoreCase',true)
            % Calculate alternans series with original beat timing
            [TWARes.lomb.at_lead(lead).series(:, :), TWARes.lomb.at_lead(lead).times] = ...
                CalcAltSeriesForLomb(ecg(:, lead), Align.fid + Align.f2s, Align.amp(:, lead), Align.st, Align.valid(:, lead));

            % Calculate Lomb-Scargle periodogram for each ST segment timepoint
            global beats
            frequencyRange = [0:beats / 2] / beats;  % 0 to 0.5 cycles/beat
            
            for timepoint = 1:size(TWARes.lomb.at_lead(lead).series, 2)
                % Detrend the alternans series for current timepoint
                detrendedSeries = detrend(TWARes.lomb.at_lead(lead).series(:, timepoint));
                
                % Compute Lomb-Scargle normalized periodogram
                [lombPSD, ~] = lomb_twa(TWARes.lomb.at_lead(lead).times', detrendedSeries, frequencyRange);
                
                % Normalize PSD to maintain energy equivalence with standard methods
                signalEnergy = sum(detrendedSeries .* detrendedSeries);
                normalizationFactor = signalEnergy / sum(lombPSD) / length(detrendedSeries);
                TWARes.lomb.psd(:, lead, timepoint) = lombPSD * normalizationFactor;
            end
        end
    
        %% Methods 2 & 3: Replacement and Difference approaches
        replacementMethods = {'replace', 'differences'};
        for methodIdx = 1:length(replacementMethods)
            methodName = replacementMethods{methodIdx};
            
            if contains(Param.MethodForEctopy, methodName, 'IgnoreCase', true)
                useDifferenceMethod = (methodIdx == 2);  % True for 'differences' method
                
                % Calculate alternans series with specified ectopy handling
                [TWARes.(methodName).at_lead(lead).series(:, :), avgEven, avgOdd] = ... 
                    CalcAltSeries(ecg(:, lead), Align.fid + Align.f2s, Align.amp(:, lead), Align.st, Align.valid(:, lead), useDifferenceMethod);
                
                % Store average even/odd beat templates if available
                if ~isempty(avgEven)
                    TWARes.(methodName).avg_even(:, lead) = avgEven;
                    TWARes.(methodName).avg_odd(:, lead) = avgOdd;
                end
            
                % Calculate power spectral density for the alternans series
                leadSeries = squeeze(TWARes.(methodName).at_lead(lead).series(:, :));
                TWARes.(methodName).psd(:, lead, :) = CalcPSD(leadSeries, useDifferenceMethod);
            end
        end

    %% Calculate final TWA significance metrics for all methods
        for methodIdx = 1:length(methodNames)
            methodName = methodNames{methodIdx};
            
            if ~isfield(TWARes, methodName)
                continue;
            end

            % Average PSD across all ST segment timepoints
            leadPSD = squeeze(TWARes.(methodName).psd(:, lead, :));
            TWARes.(methodName).avg_psd(:, lead) = CalcAvgPSD(leadPSD);
        
             % Compute TWA magnitude and statistical significance
            [TWARes.(methodName).significant(lead), TWARes.(methodName).VAlt(lead), TWARes.(methodName).Ratio(lead)] = ...
                CalcValues(TWARes.(methodName).avg_psd(:, lead));

            TWARes.(methodName).successfull = true;
        end
    
    end
end