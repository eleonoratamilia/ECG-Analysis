function Align = AlignSingleBeat(ecg, q, s, Align, CurrAvg)
% AlignSingleBeat.m
% Author: Alexander Khaustov; alexander dot khaustov at gmail dot com 
% Copyright (C) 2008 St.-Petersburg Institute of Cardiological Technics (Incart), www.incart.ru
% This software is released under the terms of the GNU General
% Public License (http://www.gnu.org/copyleft/gpl.html).
% 
% Sequential beat alignment vs current average for MMA

ind = length(Align.fidBase) + 1;

Align.fidBase(ind) = FindFidBase(ecg(:, Align.lead), q, s, Align.orientation);

[a, b, c, d] = ...
    AdjustFiducials(ecg(:, Align.lead), Align.fidBase(ind), Align.q2f, Align.f2s, CurrAvg.qs_avg, Align.st, CurrAvg.st_avg);

if (~isempty(c))
    Align.fidQRS(ind) = a;
    Align.QRScorr(ind, Align.lead) = b;
    Align.fid(ind) = c;
    Align.Tcorr(ind, Align.lead) = d;
    Align.valid(ind, 1:size(ecg, 2)) = true; %   to be done later!!!!
else
    Align.valid(ind, 1:size(ecg, 2)) = false;
end

return