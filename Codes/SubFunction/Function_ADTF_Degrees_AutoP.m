
function SU_Diff_ADTF = Function_ADTF_Degrees_AutoP(ADTF_Sub, Sub_Name, Output_file_Path)

ADTF_Sub = permute(ADTF_Sub, [2 3 1 4]);
ADTF_Sub = sum(ADTF_Sub, 4, 'omitnan') / (size(ADTF_Sub, 4)-1);
T = eye(size(ADTF_Sub, 1), size(ADTF_Sub, 2));
T = repmat(T, 1, 1, size(ADTF_Sub, 3));
ADTF_Sub(logical(T)) = 0;

ADTF_SubRep = reshape(ADTF_Sub, size(ADTF_Sub, 1), size(ADTF_Sub, 2), 500, []);
Diff_ADTF = sum(ADTF_SubRep, 2, 'omitnan');
Diff_ADTF = mean(Diff_ADTF, 3, 'omitnan');
Diff_ADTF = sum(Diff_ADTF, 1, 'omitnan');
SU_Diff_ADTF = squeeze(Diff_ADTF);

SU_Diff_ADTF_A = SU_Diff_ADTF(:);
SU_Diff_ADTF_A(SU_Diff_ADTF_A == 0) = [];
SU_Diff_ADTF_A = zscore(SU_Diff_ADTF_A);
SU_Diff_ADTF(SU_Diff_ADTF~=0) = SU_Diff_ADTF_A;

Output_file_Sub_Path = fullfile(Output_file_Path, 'Degrees');
if ~isfolder(Output_file_Sub_Path)
    mkdir(Output_file_Sub_Path);
end

Save_Mat_Path = fullfile(Output_file_Sub_Path, ['SU_Diff_', Sub_Name, '.mat']);

save(Save_Mat_Path, 'SU_Diff_ADTF');
end


function tmpDiff_ADTF = function_CalDelB_TRs(ADTF_Sub)
% 删除前n秒不稳定数据
Diff_ADTF = sum(ADTF_Sub, [1, 2], 'omitnan');
Diff_ADTF = squeeze(Diff_ADTF);
Diff_ADTF = Diff_ADTF(:);

XC = Diff_ADTF(50*100:300*100);
XM = mean(XC);
XS = std(XC);
Y3 = Diff_ADTF < (XM-4*XS);
Y4 = Diff_ADTF > (XM+4*XS);
tmpDiff_ADTF = Y3 | Y4;
% tmpDiff_ADTF = ~tmpDiff_ADTF;

end
