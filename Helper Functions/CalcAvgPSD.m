function avg_psd = CalcAvgPSD(psd)

avg_psd = mean(psd, 2); % average over all beats within the analysis window

return
