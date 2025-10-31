function [TWARes, index, end_times] = MMA_twa(ecg, freq, Def_Params)
%MMA_TWA Calculate T-Wave Alternans using Modified Moving Average method.
%   [TWARes, index, end_times] = MMA_TWA(ecg, freq, Def_Params) computes
%   T-wave alternans using the MMA algorithm on the input ECG signal.
%
%   Inputs
%   ------
%   ecg        : numeric vector
%       ECG signal data
%   freq       : positive scalar  
%       Sampling frequency in Hz
%   Def_Params : struct
%       Parameter structure containing MMA configuration fields
%
%   Outputs
%   ------ 
%   TWARes    : struct
%       TWA analysis results structure
%   index     : numeric vector
%       Start indices of TWA segments
%   end_times : numeric vector  
%       End indices of TWA segments
    
    %% Initialize global parameters
    clear global Param
    global Param beats lead
    
    Param.Metric = 'MMA';       % Modified moving average
    Param.Alignment = 'st';     % ST alignment method
    beats = Def_Params.beats;

    %% Step 1: ECG annotation and preprocessing
    ann = Annotation_Table(ecg, freq);
    [s, tend, q, qs, stend] = Interpret_Annotations(ann, ecg);
    
    % Approximate ST segment length for analysis
    stlen = ApproximateSTLen(s, stend);
    
    % Apply median filtering for noise reduction
    ecg = FilterForTWA(ecg, freq);

    %% Step 2: Configure MMA parameters from definition structure
    Param.stAdjIntv = floor(Def_Params.MMA.stAdjIntv_Multiplier * freq);
    Param.Interval = Def_Params.MMA.Interval;     % Result output frequency
    Param.corrQRS = Def_Params.MMA.corrQRS;       % QRS correlation threshold
    Param.corrT = Def_Params.MMA.corrT;           % T-wave correlation threshold  
    Param.NInapr = Def_Params.MMA.NInapr;         % Beat acceptance threshold

    %% Step 3: Perform TWA calculation using MMA method
    [TWARes] = TWAMMA_calculator(Param, ecg, q, s, stlen, freq);
    
    %% Step 4: Extract and format results
    index = [];
    end_times = [];
    
    if isfield(TWARes, 'res')
        index = s([TWARes.res.StartQRSInd]);
        end_times = s([TWARes.res.EndQRSInd]);
    end
end