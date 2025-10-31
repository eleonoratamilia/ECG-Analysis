function [CurrAvg, Align] = InitializeMMAavg(ecg, q, s, stlen)
% InitializeMMAavg.m
% Author: Alexander Khaustov; alexander dot khaustov at gmail dot com 
% Copyright (C) 2008 St.-Petersburg Institute of Cardiological Technics (Incart), www.incart.ru
% This software is released under the terms of the GNU General
% Public License (http://www.gnu.org/copyleft/gpl.html).
% 
% initialize averages for MMA

CurrAvg = [];

% first reasonable beats to start even and odd averaging
    [beatindex, Align] = FindFirstBeatsMMA(ecg, q, s, stlen);
    if beatindex == -1
        disp('TWAbyMMAOnAFile: failed to find reasonable template for even and odd beats to start averaging');
        return;
    end
    
    leads = size(ecg, 2);
    
    for i = 1:2
        for lead = 1:leads
            CurrAvg(i, lead).qs_avg = Align.qs_templ(:, lead);
            CurrAvg(i, lead).st_avg = Align.st_templ(:, lead);
        end
    end

return;