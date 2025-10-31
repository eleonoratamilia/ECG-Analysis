function arr = TWAUpdateArray(arr, beats)
%TWAUPDATEARRAY Replicate last two beats to reach target beat count for spectral analysis.
%   arr = TWAUPDATEARRAY(arr, beats) extends beat measurement arrays to the
%   required length (typically 128 beats) by replicating the last two elements
%   in an alternating pattern. Used when 90-100% of required beats are available.
%
%   Inputs
%   ------
%   arr   : numeric vector
%       Beat measurement array (QRS positions, intervals, etc.)
%   beats : positive integer
%       Target beat count for spectral TWA analysis (typically 128)
%
%   Output
%   ------
%   arr   : numeric vector
%       Extended array of length 'beats' with replicated elements
%
%   Usage Context:
%     Called when 90-100% of required beats are available to reach exact
%     beat count needed for spectral method analysis while maintaining
%     even-odd beat pattern integrity.

    currentLength = length(arr);
    
    % Extend array if shorter than required beat count
    if currentLength < beats
        % Fill odd positions with second-to-last element
        arr(currentLength + 1 : 2 : beats) = arr(currentLength - 1);
        
        % Fill even positions with last element  
        arr(currentLength + 2 : 2 : beats) = arr(currentLength);
    end
end