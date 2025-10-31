function [beatindex, Align] = FindFirstBeatsMMA(ecg, q, s, stlen)
% FindFirstBeatsMMA.m
% Author: Alexander Khaustov; alexander dot khaustov at gmail dot com 
% Copyright (C) 2008 St.-Petersburg Institute of Cardiological Technics (Incart), www.incart.ru
% This software is released under the terms of the GNU General
% Public License (http://www.gnu.org/copyleft/gpl.html).
% 
% picking first beat as initial template is dangerous - it might be very
% noisy, so we find a beat that has good correlation with majority of first
% 50 beats and therefore has acceptable noise level

    indQRS = 1;
    beats = 50; % overkill is better than agony
    
    while (indQRS + beats < length(q))
       Align = AlignBeats(ecg, beats, q(indQRS:indQRS + beats - 1), s(indQRS:indQRS + beats - 1), stlen);
       if (~isempty(Align.fid))
           break;
       end
       indQRS = indQRS + beats;
    end
   
   if (~exist('Align', 'var')) || (~isfield(Align, 'fid')) ||  (isempty(Align.fid))
       beatindex = -1;
       Align.fid=[];
   else
       ind = indQRS + Align.template - 1;
       ev_odd = mod(ind, 2) + 1;
       beatindex(ev_odd) = ind;

       for i = Align.template + 1:2:beats
           if (Align.valid(i, Align.lead))
               break;   % will break anyway, since we've got no more than 10% invalid beats
           end
       end   
       odd_ev = 3 - ev_odd;
       beatindex(odd_ev) = indQRS + i - 1;
   end
    
return