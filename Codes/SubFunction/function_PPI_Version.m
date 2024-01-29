
function [beta, bint, r, rint, stats] = function_PPI_Version(Bold_Y, Bold_X, PsySignal, TR)

[hrf, ~] = spm_hrf(TR);
PsySignal_HRF = spm_Volterra(struct('u',PsySignal,'name',{{'PsySignal'}}), hrf);

ROI_neuro = deconvwnr(Bold_X',hrf);

%%%%%%%%%%%%%%%%%%%%%% PPI analysis %%%%%%%%%%%%%%%%%%%%%%%%%
B_inter = spm_Volterra(struct('u',ROI_neuro.*PsySignal,'name',{{'RP'}}), hrf);

B_inter = B_inter(length(hrf):end);
Bold_Y = Bold_Y(length(hrf):end);
Bold_X = Bold_X(length(hrf):end);
PsySignal_HRF = PsySignal_HRF(length(hrf):end);

X = [B_inter, Bold_X, PsySignal_HRF, ones(size(PsySignal_HRF, 1), 1)];
[beta, bint, r, rint, stats] = regress(Bold_Y, X);

end


