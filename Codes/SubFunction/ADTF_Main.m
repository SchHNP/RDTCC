

Input_file_Path = '/Data/EEG';
Output_file_Path = '/ADTF';

ADTD_freqLow = 1;
ADTD_freqUpper = 45;
Channel = {'FP1' 'FP2' 'FPZ' 'F7' 'F8' 'F3' 'F4' 'FZ' 'T7' 'T8' 'C3' 'C4' 'CZ' 'P7' 'P8' 'P3' 'P4' 'PZ' 'O1' 'O2' 'OZ'};

Load_file = dir([Input_file_Path, filesep, '*.set']);
Subject_Num = size(Load_file, 1);

parfor Sub_i = 1:Subject_Num

    EEG = pop_loadset('filename', Load_file(Sub_i).name, 'filepath', Load_file(Sub_i).folder);
    [EEG, ~] = pop_select( EEG, 'channel', Channel );

    EEG = pop_resample(EEG, 100);

    fs = EEG.srate;

    Data_EEG = EEG.data;
    Data_EEG = double(Data_EEG');

    [bic, aic] = mvar_model_order1(Data_EEG', 2, 20);

    if ~isfolder(Output_file_Path)
        mkdir(Output_file_Path);
    end
    [~, Sub_Name, ~] = fileparts(Load_file(Sub_i).name);
    Sub_Name = Sub_Name(4:25);

    ADTF_Sub =  ADTF(Data_EEG, ADTD_freqLow, ADTD_freqUpper, bic, fs, 1);

    fprintf('%3d | %3d\n', Sub_i, Subject_Num)
    function_save_ADTF(ADTF_Sub, Output_file_Path, Sub_Name, ADTD_freqLow, ADTD_freqUpper);


end

function function_save_ADTF(ADTF_Sub, Output_file_Path, Sub_Name, ADTD_freqLow, ADTD_freqUpper)

subSample_out = [Output_file_Path, filesep, 'ADTF_Result'];
if ~isfolder(subSample_out)
    mkdir(subSample_out);
end

out_sub_path = fullfile(subSample_out, ['ADTF_freq_',...
    strrep(sprintf('%2d', ADTD_freqLow), ' ', '0'), '_', strrep(sprintf('%2d', ADTD_freqUpper), ' ', '0'), '_',...
    Sub_Name, '.mat']);

save(out_sub_path, 'ADTF_Sub', '-v7.3');

end


