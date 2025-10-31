function ecgout = FilterForTWA(ecgin, freq)
% FilterForTWA.m
% Author: Alexander Khaustov; alexander dot khaustov at gmail dot com 
% Copyright (C) 2008 St.-Petersburg Institute of Cardiological Technics (Incart), www.incart.ru
% This software is released under the terms of the GNU General
% Public License (http://www.gnu.org/copyleft/gpl.html).
% 
% right now only simple median filtration's used

% ecgout = ecgin;
% return;

disp('TWA: Signal filtration...');

for i = 1:size(ecgin, 2)
    ecgout(:, i) = medianfilter(ecgin(:, i), freq);
end 

disp('TWA: done');

return