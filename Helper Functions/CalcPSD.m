function psd = CalcPSD(x, calc_from_diff)
% CalcPSD.m
% Author: Alexander Khaustov; alexander dot khaustov at gmail dot com 
% Copyright (C) 2008 St.-Petersburg Institute of Cardiological Technics (Incart), www.incart.ru
% This software is released under the terms of the GNU General
% Public License (http://www.gnu.org/copyleft/gpl.html).
% 
% Calculates the Power Spectral Density of TWA series for each point
% withing ST segment
psd = zeros(floor(size(x,1)/2)+1, size(x, 2));

global Param

for timept = 1:size(x, 2)
    a = x(:, timept);
    if ~calc_from_diff
        a = detrend(a);
    end
    len = length(a);

    b = periodogram(a,hamming(len),len);
    psd(:, timept) = b * sum(a .* a) / sum(b) / length(a);    % normalization independent of periodogram estimation method
end
return
