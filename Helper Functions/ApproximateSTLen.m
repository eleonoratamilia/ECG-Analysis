function stlen = ApproximateSTLen(s, stend)
% ApproximateSTLen.m
% Author: Alexander Khaustov; alexander dot khaustov at gmail dot com 
% Copyright (C) 2008 St.-Petersburg Institute of Cardiological Technics (Incart), www.incart.ru
% This software is released under the terms of the GNU General
% Public License (http://www.gnu.org/copyleft/gpl.html).
% 
% average ST interval estimation

if (max(stend) == 0)
    rr = mean(s(2:length(s)) - s(1:length(s) - 1));
    stlen = ceil(0.5 * rr);
else
    stlen = floor(median(stend(find(stend))));
end
    
return