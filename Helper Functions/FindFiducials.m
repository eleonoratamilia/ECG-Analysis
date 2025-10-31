function Align = FindFiducials(ecg, q, s, Align)
% FindFiducials.m
% Author: Alexander Khaustov; alexander dot khaustov at gmail dot com 
% Copyright (C) 2008 St.-Petersburg Institute of Cardiological Technics (Incart), www.incart.ru
% This software is released under the terms of the GNU General
% Public License (http://www.gnu.org/copyleft/gpl.html).
% 
% finds the lead with maximum QRS amplitude, takes a beat as QRS template and finds fiducial points such that
% aligned beats have maximum cross-correlation (over Q to S or S to T interval) in the chosen
% lead against the template
%
% no attempt has been made to make this code computationally effective

global Param
q = q(1:min(length(q), length(s)));
s = s(1:min(length(q), length(s)));

leadSorted = SortLeads(ecg, q, s, Align.st);

for ilead = 1:length(leadSorted)
    Align.lead = leadSorted(ilead);

    Align.orientation = MaxFrontOrientation(ecg(:, Align.lead), q, s);

    Align.fidBase = FindFidBase(ecg(:, Align.lead), q, s, Align.orientation);


    Align.f2s = floor(median(s - Align.fidBase) + Param.stAdjIntv);   % reserve to avoid including part of QRS in ST interval computations
    Align.q2f = floor(median(Align.fidBase - q));
    Align.fid = [];
    NFailed = 0;

    Align.template = 0;
    while isempty(Align.fid)
        %    template = floor(rand * length(s)) + 1;
        Align.template = Align.template + 1;
        %disp(['template beat - ' num2str(Align.template)]);

        Align = Adjust(ecg, Align);
        if isempty(Align.fid)

            %       disp('doesn''t do...');
            NFailed = NFailed + 1;
            if (NFailed > Param.NInapr * length(s))
                         disp('what a mess..');
                break;
            end
        end
    end
    if (~isempty(Align.fid))
        for i = 1:length(Align.fid)
            for lead=1:size(ecg, 2)
                Align.amp(i, lead) = mean(ecg(Align.fidQRS(i) - Align.q2f : Align.fidQRS(i) + Align.f2s, lead));  % align to zero mean on QS segment
            end
        end
        Align = CheckLeadCorrelations(ecg, Align);
        disp(['FindFiducials: Alignment in lead ' num2str(Align.lead)]);
        break;
    end
    disp(['FindFiducials: lead ' num2str(Align.lead) ' is too noisy or ectopic..']);
end
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function leads = SortLeads(ecg, q, s, st)

NLead = size(ecg, 2);

global Param
%   QRS or T amplitudes in each lead
ptp = zeros(min(length(q), 8),NLead);
for i = 1:min(length(q), 8)
    for l = 1:NLead
        if strcmp(Param.Alignment, 'st')
            ptp(i, l) = max(ecg(s(i):s(i) + st, l)) - min(ecg(s(i):s(i) + st, l));
        else
            ptp(i, l) = max(ecg(q(i):s(i), l)) - min(ecg(q(i):s(i), l));
        end
    end
end
ptpmed = zeros(1,NLead);
for i = 1:NLead
    ptpmed(i) = median(ptp(:, i));
end

[sorted indices] = sort(ptpmed);
leads = indices(length(indices):-1:1);
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Align = Adjust(ecg, Align)

global Param

qs_templ = ecg((Align.fidBase(Align.template) - Align.q2f):(Align.fidBase(Align.template) + Align.f2s), Align.lead);
st_templ = ecg((Align.fidBase(Align.template) + Align.f2s):(Align.fidBase(Align.template) + Align.f2s + Align.st), Align.lead);

[a b c d] = AdjustFiducials(ecg, Align.fidBase, Align.q2f, Align.f2s, qs_templ, Align.st, st_templ);

if ~isempty(c)
    Align.fidQRS = a;
    Align.QRScorr(:, Align.lead) = b;
    Align.fid = c;
    Align.Tcorr(:, Align.lead) = d;
    amp = mean(qs_templ);
    Align.qs_templ(:, Align.lead) = qs_templ - amp;
    Align.st_templ(:, Align.lead) = st_templ - amp;
    
    for i = 1:length(Align.fid)
        Align.valid(i, Align.lead) = (Align.QRScorr(i, Align.lead) >= Param.corrQRS) && (Align.Tcorr(i, Align.lead) >= Param.corrT);
    end
    Align.validleads(Align.lead) = true;
end

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function orientation = MaxFrontOrientation(ecg, q, s)

NPos = 0;

for i = 1:min(20, length(s))
    [ma maI] = max(ecg(q(i):s(i)));
    [mi miI] = min(ecg(q(i):s(i)));

    NPos = NPos + (maI > miI);
end

orientation = (NPos > 10);

return
