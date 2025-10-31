function flatline = check_flatline(sig)
% CHECK_FLATLINE checks for flatline segments by partitioning the signal 
% into 10-sample frames and testing if max == min in any frame.
%
%   Inputs
%   ------
%   sig : numeric vector
%       Input signal (expected to be partitioned/compatible with buffering).
%
%   Output
%   ------
%   flatline : logical
%       True if any 10-sample segment has constant value (flatline detected).

    % Partition signal into 10-sample frames
    ba = buffer(sig, 10);
    
    % Check for flatlines (max == min in any segment)
    flatline = any(max(ba) == min(ba));
end