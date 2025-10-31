function fidBase = FindFidBase(ecg, q, s, orientation)
% FindFidBase.m
% Author: Alexander Khaustov; alexander dot khaustov at gmail dot com 
% Copyright (C) 2008 St.-Petersburg Institute of Cardiological Technics (Incart), www.incart.ru
% This software is released under the terms of the GNU General
% Public License (http://www.gnu.org/copyleft/gpl.html).
% 
% basic fiducial point - middle of the maximum front between Q and S; to be
% improved by adjustment later

fidBase = zeros(1,length(s));
for i = 1:length(s)
    [startI endI] = FindFront(ecg(q(i):s(i)), orientation);
    startI = startI + q(i) - 1;
    endI = endI + q(i) - 1;
%    [m fidBase(i)] = max(abs(alt_series(ecg(startI:endI))));     % quite unstable, changed that
%    fidBase(i) = fidBase(i) + startI - 1;
    half = 0.5 * abs(ecg(endI) - ecg(startI));
    for j = startI:endI        
        if abs(ecg(j) - ecg(startI)) > half
            % closer to half between j and j - 1
            if (half - abs(ecg(j - 1) - ecg(startI)) > abs(ecg(j) - ecg(startI)) - half)
                fidBase(i) = j;
            else
                fidBase(i) = j - 1;
            end
            break;
        end
    end
    
    if (fidBase(i) == 0)
        fidBase(i) = startI;
    end
end

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [startI, endI] = FindFront(ecg, orientation)
% front with maximum amplitude and desired orientation
[ma maI] = max(ecg);
[mi miI] = min(ecg);

if (maI == miI)
    startI = 1;
    endI = length(ecg);
    return;
end

if (maI > miI) == orientation
    startI = min(miI, maI);
    endI = max(miI, maI);
elseif (orientation == false)
    [pre preI] = max(ecg(1:miI));
    [post postI] = min(ecg(maI:length(ecg)));
    if (pre - mi > ma - post)
        startI = preI;
        endI = miI;
    else
        startI = maI;
        endI = postI + maI - 1;
    end
else
    [pre preI] = min(ecg(1:maI));
    [post postI] = max(ecg(miI:length(ecg)));
    if (ma - pre > post - mi)
        startI = preI;
        endI = maI;
    else
        startI = miI;
        endI = postI + miI - 1;
    end
end

if (startI == endI)
    startI = min(miI, maI);
    endI = max(miI, maI);
end
return
