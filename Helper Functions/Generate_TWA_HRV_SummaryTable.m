function final_table = Generate_TWA_HRV_SummaryTable(ecg_count, T_spectral_all, T_MMA_all, RR_all)
%GENERATESUMMARYTABLE Create final summary table with TWA and HRV metrics
    
    final_table = table('Size', [1, 9], 'VariableTypes', repmat({'double'}, 1, 9), ...
        'VariableNames', {'Ecg_count', 'Spec_Max', 'MMA_Max', 'Spec_Median', 'MMA_Median', ...
                          'HR', 'HRV_SDNN', 'HRV_RMSSD', 'HRV_rrHRV'});
    
    final_table.Ecg_count(1) = ecg_count;
    
    % HR & HRV
    if ~isempty(RR_all)
        RR = HRV.RRfilter(RR_all, 0.15);
        final_table.HR = HRV.HR(RR);
        final_table.HRV_SDNN = HRV.SDNN(RR);
        final_table.HRV_RMSSD = HRV.RMSSD(RR);
        final_table.HRV_rrHRV = HRV.rrHRV(RR);
    else
        final_table.HR = NaN;
        final_table.HRV_SDNN = NaN;
        final_table.HRV_RMSSD = NaN;
        final_table.HRV_rrHRV = NaN;
    end
    
    if isempty(T_spectral_all) || ~any(T_spectral_all.Ratio >= 3)
        final_table.Spec_Max = NaN;
        final_table.Spec_Median = NaN;
    else
        spectral = T_spectral_all.Value;
        spectral = spectral(T_spectral_all.Ratio >= 3);
        final_table.Spec_Max = nanmax(spectral);
        final_table.Spec_Median = nanmedian(spectral);
    end
    
    % TWA MMA
    if isempty(T_MMA_all)
        final_table.MMA_Max = NaN;
        final_table.MMA_Median = NaN;
    else
        mma = T_MMA_all.Value;
        mma = mma(mma ~= 0);
        final_table.MMA_Max = nanmax(mma);
        final_table.MMA_Median = nanmedian(mma); 
    end
end