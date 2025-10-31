function [s, tend, q, qs, stend] = Interpret_Annotations(ann, ecg)
%INTERPRET_ANNOTATIONS Extract QRS complex and ST segment features from annotations.
%   [s, tend, q, qs, stend] = INTERPRET_ANNOTATIONS(ann, ecg) processes WFDB-style
%   annotations to extract Q-wave onset, S-wave position, T-wave end, and
%   derived intervals.
%
%   Inputs
%   ------
%   ann : struct array
%       Annotation structure with fields: time, anntyp, num, chan
%   ecg : numeric matrix
%       ECG signal data (used for bounds checking)
%
%   Outputs
%   ------
%   s     : numeric vector - S-wave positions (QRS offset)
%   tend  : numeric vector - T-wave end positions  
%   q     : numeric vector - Q-wave onset positions
%   qs    : numeric vector - QRS durations (s - q)
%   stend : numeric vector - ST segment durations (tend - s)

% Note : Original file mentioned the contribution of Shamim

    %% Initialize output variables
    s = []; tend = []; q = []; qs = []; stend = [];
    qlatest = -1;  % Tracks most recent Q-wave onset
    
    %% Process annotations sequentially
    i = 1;
    while i <= length(ann)
        % Q-wave onset annotation (type 39, num 1)
        if (ann(i).num == 1 && ann(i).anntyp == 39 && ann(i).chan == 0)
            if ann(i).time <= size(ecg, 1)
                qlatest = ann(i).time;
            end
            
        % S-wave annotation (type 40, num 1) - marks QRS end
        elseif (ann(i).num == 1 && ann(i).anntyp == 40 && ann(i).chan == 0)
            if ann(i).time <= size(ecg, 1)
                % Valid Q-S pair found: store S position and calculate QS interval
                if (qlatest ~= -1 && (isempty(s) || qlatest > s(end)))
                    s(end + 1) = ann(i).time;
                    q(end + 1) = qlatest;
                    qs(end + 1) = s(end) - q(end);
                end

                % Calculate ST interval for previous S-wave if T-end available
                if length(s) > 1
                    if isempty(tend) || tend(end) < s(end - 1)
                        stend(end + 1) = 0;  % No valid T-end for this S-wave
                    else
                        stend(end + 1) = tend(end) - s(end - 1);
                    end
                end
            end
            
        % T-wave end annotation (type 40, num 2)
        elseif (ann(i).num == 2 && ann(i).anntyp == 40 && ann(i).chan == 0)
            if ann(i).time <= size(ecg, 1)
                tend(end + 1) = ann(i).time;
            end
        end
        
        i = i + 1;
    end

    %% Process final ST segment if T-end available for last S-wave
    if ~isempty(tend) && tend(end) > s(end)
        stend(end + 1) = tend(end) - s(end);
    end

    %% Remove last QRS complex to avoid boundary issues
    s = s(1:end - 1);
    q = q(1:length(s));
end