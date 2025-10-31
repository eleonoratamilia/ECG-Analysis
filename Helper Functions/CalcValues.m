%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [significant, VAlt, Ratio] = CalcValues(avg_psd)
% CalcValues.m
% Author: Alexander Khaustov; alexander dot khaustov at gmail dot com 
% Copyright (C) 2008 St.-Petersburg Institute of Cardiological Technics (Incart), www.incart.ru
% This software is released under the terms of the GNU General
% Public License (http://www.gnu.org/copyleft/gpl.html).
% 
% CalcValues calculates the alternans values for Spectral Method
% INPUT:
%    avg_psd (num_of_beats/2+1 x num_of_leads)  is the psd averaged across the entire S-Tend segment
% OUTPUT:
%    significant  (1 x num_of_leads) 1 if the alternans value is statistically significant against the Param.RatioThreshold in the lead and 0 otherwise
%    VAlt   (1 x num_of_leads) spectral amplitude @ 0.5 cycles-per-beat -
%      'voltage of alternans'
%    Ratio  (1 x num_of_leads) 'alternans ratio'

global Param

NoiseStart = floor(2 * 0.4 * (length(avg_psd) - 1) + 1);
NoiseEnd = floor(2 * 0.46 * (length(avg_psd) - 1) + 1);


% Calculate the mean and std of the noise
nmean = mean(avg_psd(NoiseStart : NoiseEnd - 1));
nstd = std(avg_psd(NoiseStart : NoiseEnd - 1));
% spectral amplitude @ 0.5 cycles-per-beat
v = avg_psd(length(avg_psd));
if ((v - nmean) > 0)
    VAlt = sqrt(v - nmean);
else
    VAlt=-1;
end
    Ratio = (v - nmean) / nstd;

significant = 1;

return
