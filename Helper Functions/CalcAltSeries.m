function [alt_s, avg_even, avg_odd] = CalcAltSeries(ecg, fid, amp, STLen, valid, calc_from_diff)
% CalcAltSeries.m
% Author: Alexander Khaustov; alexander dot khaustov at gmail dot com 
% Copyright (C) 2008 St.-Petersburg Institute of Cardiological Technics (Incart), www.incart.ru
% This software is released under the terms of the GNU General
% Public License (http://www.gnu.org/copyleft/gpl.html).
% 
% calculates the alternans series for different points within st segment. 
% If a given beat is marked as invalid it's replaced by the average
% even/odd beat
% this is the routine used for Param.MethodForEctopy = 'replace' and
% Param.MethodForEctopy = 'differences'
%
% normally the series is formed from sequential amplitudes of equally timed
% points within ST segment, but if the parameter calc_from_diff is
% set - series is calculated as differences between amplitudes in adjacent
% beats which performs natural detrending of the data that is later subject
% to power spectral estimation and supresses low frequencies
%
% INPUT:
%       see AlignBeats.m
% OUTPUT: 
%       alt_s       (num_of_beats x num_of_leads x num_of_timepoints): TWA series                 
%       avg_even (num_of_timepoints x num_of_leads)
%       avg_odd  (num_of_timepoints x num_of_leads)

    avg_even = [];
    avg_odd = [];
    invalid_exist = ~isempty(find(valid==0));
    if (invalid_exist)
        vinds = find(valid);
        odd = vinds(find(mod(vinds(:), 2)));  %   odd indices of valid QRS
        even = vinds(find((1 - mod(vinds(:), 2))));   %   even indices of valid QRS
    end    
    alt_s = zeros(length(fid) - 1, STLen);

    for timept = 1:STLen           
        %   whenever there are invalid complexes compute odd and even
        %   averages to replace those
        if (invalid_exist)
            avg_odd(timept) = mean(ecg(fid(odd) + timept) - amp(odd));                
            avg_even(timept) = mean(ecg(fid(even) + timept) - amp(even));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for i = 1:(length(fid) - 1)
            
            if (valid(i)) 
                first = ecg(fid(i) + timept) - amp(i);
            elseif mod(i, 2)
                first = avg_odd(timept);
            else
                first = avg_even(timept);
            end

            if (~calc_from_diff)
                alt_s(i, timept) = first;
            else
                if (valid(i + 1)) 
                    second = ecg(fid(i + 1) + timept) - amp(i + 1);
                elseif mod(i + 1, 2)
                    second = avg_odd(timept);
                else
                    second = avg_even(timept);
                end
                alt_s(i, timept) = first - second;
            end
            
        end
    end  

return
