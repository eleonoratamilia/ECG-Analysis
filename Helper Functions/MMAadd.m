function [avg_new, invalid] = MMAadd(avg_old, v)
% MMAadd.m
% Author: Alexander Khaustov; alexander dot khaustov at gmail dot com 
% Copyright (C) 2008 St.-Petersburg Institute of Cardiological Technics (Incart), www.incart.ru
% This software is released under the terms of the GNU General
% Public License (http://www.gnu.org/copyleft/gpl.html).
% 
% Threshold averaging for MMA

invalid = 0;
for i = 1:length(v)
    df = (v(i) - avg_old(i)) / 8;
    if (abs(df) >= 32) % must be mV?
        avg_new(i) = avg_old(i) + sign(df) * 32;
        invalid = invalid + 1;
    elseif (abs(df) >= 1)
        avg_new(i) = avg_old(i) + df;
    else
        avg_new(i) = avg_old(i) + sign(df);
    end    
end

return