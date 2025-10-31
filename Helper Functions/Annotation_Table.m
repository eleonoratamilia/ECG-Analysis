function ann = Annotation_Table(ecg, freq)
%ANNOTATIONTABLE Convert ECG fiducial points to WFDB-style annotations.
%   ann = ANNOTATIONTABLE(ecg, fs) runs Annotate_ECG_Multi to obtain the
%   Fiducial Point Table (FPT), cleans edge cases, and returns a struct array
%   of annotations with fields: time, anntyp, subtyp, chan, num, aux.
%
%   Inputs
%   ------
%   ecg : numeric vector
%       Single-channel ECG signal.
%   freq  : positive scalar
%       Sampling frequency in Hz (passed through to Annotate_ECG_Multi).
%
%   Output
%   ------
%   ann : struct (Nx1)
%       Chronologically sorted annotations:
%         - QRS onset  -> anntyp = 39, num = 1  (FPT col 4)
%         - QRS offset -> anntyp = 40, num = 1  (FPT col 8)
%         - T offset   -> anntyp = 40, num = 2  (FPT col 12)

    %% Get Fiducial Point Table
    [~, FPT_Cell] = Annotate_ECG_Multi(ecg, freq); 
    FPT = FPT_Cell{1, 1}; 
    n_sample = length(ecg);

    %% Clean edges / invalid beats
    % Remove last beat if incomplete (T-wave end beyond signal length)
    if isempty(FPT(end, 12)) || FPT(end, 12) > n_sample
        FPT(end, :) = [];
    end
    
    % Remove first beat if P-wave is missing or invalid
    if isempty(FPT(1, 1)) || FPT(1, 1) < 2
        FPT(1, :) = [];
    end

    %% Build annotation list
    ann = struct('time', 0, 'anntyp', 0, 'subtyp', 0, 'chan', 0, 'num', 0, 'aux', '');
    k = 1;
    
    for i = 1:size(FPT, 1)
        % QRS onset annotation (FPT column 4)
        if FPT(i, 4) && FPT(i, 4) < n_sample
            ann(k).time   = FPT(i, 4);
            ann(k).anntyp = 39;  % QRS onset marker
            ann(k).subtyp = 0;
            ann(k).chan   = 0;
            ann(k).num    = 1;
            ann(k).aux    = [];
            k = k + 1;
        end
        
        % QRS offset annotation (FPT column 8)  
        if FPT(i, 8) && FPT(i, 8) < n_sample
            ann(k).time   = FPT(i, 8);
            ann(k).anntyp = 40;  % QRS/T-wave marker
            ann(k).subtyp = 0;
            ann(k).chan   = 0;
            ann(k).num    = 1;
            ann(k).aux    = [];
            k = k + 1;
        end
        
        % T-wave offset annotation (FPT column 12)
        if FPT(i, 12) && FPT(i, 12) < n_sample
            ann(k).time   = FPT(i, 12);
            ann(k).anntyp = 40;  % QRS/T-wave marker
            ann(k).subtyp = 0;
            ann(k).chan   = 0;
            ann(k).num    = 2;
            ann(k).aux    = [];
            k = k + 1;
        end
    end
    
    ann = ann(:);
end