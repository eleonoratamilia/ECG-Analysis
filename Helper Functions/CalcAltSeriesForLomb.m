function [df, times] = CalcAltSeriesForLomb(ecg, fid, amp, STLen, valid)
% CalcAltSeriesForLomb.m
% Author: Alexander Khaustov; alexander dot khaustov at gmail dot com 
% Copyright (C) 2008 St.-Petersburg Institute of Cardiological Technics (Incart), www.incart.ru
% This software is released under the terms of the GNU General
% Public License (http://www.gnu.org/copyleft/gpl.html).
% 
% calculates the alternans series for different points within st segment. 
% If a given beat is marked as invalid it is excluded and corresponding
% 'times' accounts for it
% this is the routine used for Param.MethodForEctopy = 'lomb'
%
% INPUT:
%       see AlignBeats.m
% OUTPUT: 
%       df       (number_of_valid_beats x num_of_leads x num_of_timepoints): TWA series 
%                 array (difference of even and odd beats within the nalysis window)
%       times   contains respective indices of beats where df is taken;
%           could be used for lomb


ind = 1; 
for i = 1:(length(fid) - 1)
    if (valid(i))
        times(ind) = i;
        for timept = 1:STLen
            df(ind, timept) = ecg(fid(i) + timept) - amp(i);
        end
        ind = ind + 1;
    end
end

return