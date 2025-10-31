function [TWARes] = TWAMMA_calculator(Param, ecg, q, s, stlen, freq)
%TWAMMA_CALCULATOR Compute T-Wave Alternans using Modified Moving Average method.
%   TWARes = TWAMMA_CALCULATOR(Param, ecg, q, s, stlen, freq) calculates
%   T-wave alternans by maintaining separate moving averages for even and odd
%   beats and comparing ST segment differences.
%
%   Inputs
%   ------
%   Param : struct
%       Algorithm parameters including Interval, NInapr, etc.
%   ecg   : numeric matrix
%       ECG signal data (samples x leads)
%   q     : numeric vector
%       Q-wave onset positions
%   s     : numeric vector
%       S-wave positions
%   stlen : numeric scalar
%       ST segment length in samples
%   freq  : positive scalar
%       Sampling frequency in Hz
%
%   Output
%   ------
%   TWARes : struct
%       TWA results structure containing VAlt, trends, and segment information

    %% Initialize global variables and parameters
    clear global TWARes Align CurrAvg
    global Align CurrAvg

    intervalLength = Param.Interval * freq;  % Convert interval to samples
    
    %% Initialize result structure
    TWARes.successfull = false;
    TWARes.max_lead = 1;
    TWARes.max_index = 1;

    %% Initialize moving averages and alignment
    [CurrAvg, Align] = InitializeMMAavg(ecg, q, s, stlen);
    if isempty(CurrAvg)
        return;  % Initialization failed
    end

    %% Process beats sequentially
    nextInterval = intervalLength;
    startQRSInd = 1;
    leads = size(ecg, 2);
    NInvalid = zeros(1, leads);

    for i = 1:length(q)
        even_odd = mod(i, 2) + 1;  % 1 for odd beats, 2 for even beats
            
        % Align current beat if needed
        if i > length(Align.fid)
            Align = AlignSingleBeat(ecg, q(i), s(i), Align, CurrAvg(even_odd, Align.lead));
        end
        
        %% Process each lead for current beat
        for lead = 1:leads
            if Align.valid(i, lead)
                % Calculate amplitude for baseline correction
                Align.amp(i, lead) = mean(ecg(Align.fidQRS(i) - Align.q2f : Align.fidQRS(i) + Align.f2s, lead));
    
                % Update QRS segment moving average
                qs_segment = ecg(Align.fidQRS(i) - Align.q2f : Align.fidQRS(i) + Align.f2s, lead) - Align.amp(i, lead);
                [CurrAvg(even_odd, lead).qs_avg, invalidQS] = MMAadd(CurrAvg(even_odd, lead).qs_avg, qs_segment);

                % Update ST segment moving average
                st_segment = ecg(Align.fid(i) + Align.f2s : Align.fid(i) + Align.f2s + Align.st, lead) - Align.amp(i, lead);
                [CurrAvg(even_odd, lead).st_avg, invalidST] = MMAadd(CurrAvg(even_odd, lead).st_avg, st_segment);                
                
                % Track invalid segments based on noise threshold
                qs_threshold = Param.NInapr * (Align.q2f + Align.f2s + 1);
                st_threshold = Param.NInapr * Align.st;
                if invalidQS > qs_threshold || invalidST > st_threshold
                    NInvalid(lead) = NInvalid(lead) + 1;
                end
            else
                NInvalid(lead) = NInvalid(lead) + 1;
            end
        end
        
        %% Process interval results when reaching interval boundary
        if i >= length(q) || q(i + 1) > nextInterval
            intervalNum = nextInterval / intervalLength;  % Fixed variable name
            
            % Determine valid leads based on invalid beat count
            for lead = 1:leads
                invalidThreshold = Param.NInapr * (i - 1 - startQRSInd);  
                TWARes.res(intervalNum).validleads(lead) = (NInvalid(lead) < invalidThreshold);
            end
                                    
            % Store interval metadata
            TWARes.res(intervalNum).StartQRSInd = startQRSInd;  
            TWARes.res(intervalNum).EndQRSInd = i;
            
            % Calculate TWA results for this interval
            TWARes = AddMMARes(intervalNum, TWARes);
          
            % Reset for next interval
            nextInterval = nextInterval + intervalLength;
            startQRSInd = i + 1;  % Fixed variable name
            NInvalid = zeros(1, leads);
        end
        
        %% Update global TWA maximum if successful detection
        if TWARes.successfull
            [maxVals, maxIndices] = max(TWARes.VAltTrend);
            [TWARes.maxVAlt, TWARes.max_lead] = max(maxVals);
            TWARes.max_index = maxIndices(TWARes.max_lead);
        end
    end
end

function TWARes = AddMMARes(intervalNum, TWARes)
%ADDMARES Calculate TWA values for a given interval from moving averages.
%   Computes the difference between even and odd beat ST segments and
%   identifies the maximum alternans magnitude.

    global Align CurrAvg
    leads = size(CurrAvg, 2);

    TWARes.res(intervalNum).Avg = CurrAvg;

    for lead = 1:leads
        if TWARes.res(intervalNum).validleads(lead)
            % Calculate difference between even and odd ST segment averages
            st_diff = TWARes.res(intervalNum).Avg(1, lead).st_avg - ...
                      TWARes.res(intervalNum).Avg(2, lead).st_avg;
                  
            [maxAlt, TWARes.res(intervalNum).VAltPt(lead)] = max(abs(st_diff));
            
            % Store alternans magnitude (mean difference across ST segment)
            TWARes.res(intervalNum).VAlt(lead) = abs(mean(st_diff)); % ma
            TWARes.successfull = true;
        else
            TWARes.res(intervalNum).VAlt(lead) = 0;
        end
    end

    % Update TWA trend across intervals
    TWARes.VAltTrend(intervalNum, :) = TWARes.res(intervalNum).VAlt;
end