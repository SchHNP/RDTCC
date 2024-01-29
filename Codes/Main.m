
%%
Data_EEG = rand(19, 125000);
[bic, aic] = mvar_model_order(Data_EEG, 2, 20);
fs = 250;
ADTF_Sub =  ADTF(Data_EEG', 8, 13, bic, fs, 1);


%%
Function_ADTF_Degrees_AutoP(ADTF_Sub, 'Subname', 'OutputfilePath');


%%
Bold_Y = rand(250, 1);
Bold_X = rand(250, 1);
PsySignal = rand(250, 1);
beta = function_PPI_Version(Bold_Y, Bold_X, PsySignal, 2);


%%
X = rand(30, 50);
Y = rand(30, 40);
[R, T, p, df] = dcor_uc(X, Y);

