function Align = AlignBeats(ecg, beats, q, s, STLen)
% AlignBeats.m
% Author: Alexander Khaustov; alexander dot khaustov at gmail dot com 
% Copyright (C) 2008 St.-Petersburg Institute of Cardiological Technics (Incart), www.incart.ru
% This software is released under the terms of the GNU General
% Public License (http://www.gnu.org/copyleft/gpl.html).
% 
% Interprets results of FindFiducials and prints diagnostic messages
% INPUT:
% ecg   ecg signal (could be multiple channels)
% beats number of analysis beats
% q     location of q marks
% s     location of s marks
% STlen ST interval segment length
%
% OUTPUT:FindFiducials
% Align: structure containing the fields:
%       fid      (1 x beats): empty when unsuccessful, fiducial
%          points that maximize cross-correlation on ST interval
%       fidQRS   (1 x beats): empty when unsuccessful, fiducial
%          points that maximize cross-correlation on QS interval
%       amp      (beats x leads): amplitudes to subtract to align beats vertically, leads = size(ecg, 2)
%       q2f      average estimate of interval between Q and fiducial points
%       f2s      average estimate of interval between fiducial points and S
%       st       equal to STLen
%       valid    (beats x leads): 1 if the beat correlation with template is acceptable both on QS and ST 
%           intervals in the specified lead, 0 otherwise

global Param;

    Align.fid = [];
    Align.valid = [];
    Align.st = STLen;

    Align = FindFiducials(ecg, q, s, Align);
    if (isempty(Align.fid))
        disp('no lead is suitable for alignment, finishing..');
        return;
    end
    
    for lead = 1:size(ecg, 2)
        NInvalid = length(find(1 - Align.valid(:, lead)));
        if NInvalid
            fprintf('AlignBeats: %d invalid beats in lead %d\n', NInvalid, lead);
        end
    end

    if (length(Align.fid) > beats)
        Align.fid = Align.fid(1:beats);
    end
        
return
