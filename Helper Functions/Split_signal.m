function [segmentationTable] = Split_signal(good_quality, freq, Def_Params)
%SPLIT_SIGNAL Segment good-quality ECG signal into analysis windows.
%   segmentationTable = SPLIT_SIGNAL(good_quality, freq, Def_Params) divides
%   continuous good-quality segments into analysis windows based on duration
%   constraints for TWA analysis.
%
%   Inputs
%   ------
%   good_quality : logical vector
%       Binary mask indicating good quality samples (1 = good, 0 = poor)
%   freq         : positive scalar
%       Sampling frequency in Hz
%   Def_Params   : struct
%       Parameter structure containing segmentation constraints:
%         - minGoodLength: Minimum segment duration in seconds (default: 120)
%         - splitLength:   Target segment duration in seconds (default: 240)  
%         - maxGoodLength: Maximum segment duration in seconds (default: 420)
%
%   Output
%   ------
%   segmentationTable : table
%       Table with columns: PartNumber, StartSample, EndSample

    %% Extract segmentation parameters
    minGoodLength = Def_Params.Spectral.minGoodLength;  % Minimum duration: 2 minutes
    splitLength   = Def_Params.Spectral.splitLength;    % Target split: 4 minutes
    maxGoodLength = Def_Params.Spectral.maxGoodLength;  % Maximum duration: 7 minutes
    
    %% Convert time constraints to samples
    minGoodSamples = minGoodLength * freq;
    splitSamples   = splitLength * freq;
    maxGoodSamples = maxGoodLength * freq;
    
    %% Identify continuous good-quality segments
    changes = diff([0; good_quality; 0]);  % Pad to detect edge transitions
    
    startIndices = find(changes == 1);     % Start of good segments
    endIndices   = find(changes == -1) - 1; % End of good segments
    
    %% Initialize segmentation table
    segmentationTable = table([], [], [], ...
        'VariableNames', {'PartNumber', 'StartSample', 'EndSample'});
    
    partNumber = 1;

    %% Process each good-quality segment
    for i = 1:length(startIndices)
        segmentStart = startIndices(i);
        segmentEnd   = endIndices(i);
        segmentLength = segmentEnd - segmentStart + 1;
        
        % Skip segments shorter than minimum requirement
        if segmentLength <= minGoodSamples
            continue;
        end
        
        %% Split segment according to duration rules
        remainingSamples = segmentLength;
        currentStart = segmentStart;
        
        while remainingSamples > minGoodSamples
            % Determine optimal segment end position
            if remainingSamples > maxGoodSamples
                % Case 1: Long segment -> split at target length
                currentEnd = currentStart + splitSamples - 1;
                partDuration = splitSamples;
                
            elseif remainingSamples >= splitSamples
                % Case 2: Medium segment -> use full remaining length
                currentEnd = segmentEnd;
                partDuration = remainingSamples;
                
            else
                % Case 3: Short segment -> use if above minimum threshold
                currentEnd = segmentEnd;
                partDuration = remainingSamples;
            end
            
            % Add segment to output table
            segmentationTable = [segmentationTable; {partNumber, currentStart, currentEnd}];
            
            % Update for next iteration
            currentStart = currentEnd + 1;
            remainingSamples = remainingSamples - partDuration;
            partNumber = partNumber + 1;
            
        end
    end
end