

%%
[bic, aic] = mvar_model_order(Data_EEG', 2, 20);
ADTF_Sub =  ADTF(Data_EEG, ADTD_freqLow, ADTD_freqUpper, bic, fs, 1);


%%
Function_ADTF_Degrees_AutoP(ADTF_Sub(:, :, :, ADTD_freq(ADTD_freq_i, 1):ADTD_freq(ADTD_freq_i, 2)), Sub_Name_S, Output_file_PathSub);


%%
rand()
Bold_Y = ROI_GignalSub(ROI_GignalSub_i, :);
Bold_X = ROI_GignalSub(ROI_GignalSub_j, :);
PsySignal = Psy_Sub;
beta = function_PPI_Version(Bold_Y, Bold_X, PsySignal, 2);
PPI_MatSubB1(ROI_GignalSub_i, ROI_GignalSub_j) = beta(1);


%%

Bold_Y = ROI_GignalSub(ROI_GignalSub_i, :);
Bold_X = ROI_GignalSub(ROI_GignalSub_j, :);
PsySignal = Psy_Sub;
beta = function_PPI_Version(Bold_Y, Bold_X, PsySignal, 2);
PPI_MatSubB1(ROI_GignalSub_i, ROI_GignalSub_j) = beta(1);
