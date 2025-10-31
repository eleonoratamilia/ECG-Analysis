function [fidQRS, QRScorr, fidST, Tcorr] = AdjustFiducials(ecg, fidBase, q2f, f2s, qs_templ, st, st_templ)
% AdjustFiducials.m
% Author: Alexander Khaustov; alexander dot khaustov at gmail dot com 
% Copyright (C) 2008 St.-Petersburg Institute of Cardiological Technics (Incart), www.incart.ru
% This software is released under the terms of the GNU General
% Public License (http://www.gnu.org/copyleft/gpl.html).
% 
% Takes fidBase as basic fiducial point and searches for fidQRS and fidST
% - fiducial points that maximize cross-correlation of ecg intervals
% (determined with the help of q2f, f2s and st)
% located around these points with qs_templ and st_templ
% The routine also outputs respective maximum correlations

global Param
NInapr = 0;

qoff = - q2f - Param.stAdjIntv - 1;
soff = f2s - Param.stAdjIntv - 1;
toff = f2s + st - Param.stAdjIntv - 1;

warning off MATLAB:divideByZero     % to supress corrcoef warning when data is constant, in that case corrcoef is NaN and max index is any

for i = 1:length(fidBase)
        cc = zeros(1, 2 * Param.stAdjIntv + 1);
        for j = 1:(2 * Param.stAdjIntv + 1)
            qs_ecg = ecg((fidBase(i) + j + qoff):(fidBase(i) + j + soff));
            a = corrcoef(qs_templ, qs_ecg);
            cc(j) = a(1, 2);
            
            if strcmp(Param.Alignment, 'st')
                st_ecg = ecg((fidBase(i) + j + soff):(fidBase(i) + j + toff));
                b = corrcoef(st_templ, st_ecg);
                ccc(j) = b(1, 2);
            end
        end
        [m ind] = max(cc);
        indQRS = ind;
        if (strcmp(Param.Alignment, 'st'))
            [mm ind] = max(ccc);
        end

        if (m < Param.corrQRS || (strcmp(Param.Alignment, 'st') && mm < Param.corrT))
            NInapr = NInapr + 1;
        end

        if (NInapr > Param.NInapr * length(fidBase)) % previous one was 0.1
            fidST = [];
            fidQRS = [];
            QRScorr = [];
            Tcorr = [];
            break;
        else
            fidST(i) = fidBase(i) + ind - Param.stAdjIntv - 1;
            fidQRS(i) = fidBase(i) + indQRS - Param.stAdjIntv - 1;
            QRScorr(i) = m;
            if (strcmp(Param.Alignment, 'st'))
                Tcorr(i) = mm;
            end
        end
end

return