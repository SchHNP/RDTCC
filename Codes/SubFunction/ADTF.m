function [gamma2_set] = ADTF(ts,low_freq,high_freq,p,fs,mm)
% ADTF - perform Adaptive Directed Transfer Function analysis among multi-channel time series.
%
% Usage: gamma2_set = ADTF(ts,low_freq,high_freq,p,fs)
%
% Input: ts - the time series where each column is the temporal data from a
%             single channel
%        low_freq - the lowest frequency to perform DTF analysis
%        high_freq - the highest frequency to perform DTF analysis
%        p - the order of the MVAAR model
%        fs - the sampling frequency
%
% Output: gamma2_set - the DTF values from the points in the time series
%
% Description: This function performs DTF analysis on a time series using
%              an adaptive MVAR model (MVAAR). The MVAAR model generates an
%              updated coefficient matrix for each time point which is then
%              used in the DTF calculations. The output is in the form
%              gamma2_set(a,b,c,d) where a = the time point, b = the sink
%              channel, c = the source channel, d = the index for the
%              frequency value.
%
% Program Author: Christopher Wilke, University of Minnesota, USA
%
% User feedback welcome: e-mail: econnect@umn.edu
%

% License
% ==============================================================
% This program is part of the eConnectome.
%
% Copyright (C) 2010 Regents of the University of Minnesota. All rights reserved.
% Correspondence: binhe@umn.edu
% Web: econnectome.umn.edu
%
% This program is free software for academic research: you can redistribute it and/or modify
% it for non-commercial uses, under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see http://www.gnu.org/copyleft/gpl.html.
%
% This program is for research purposes only. This program
% CAN NOT be used for commercial purposes. This program
% SHOULD NOT be used for medical purposes. The authors
% WILL NOT be responsible for using the program in medical
% conditions.
% ==========================================

% Revision Logs
% ==========================================
%
% Yakang Dai, 19-June-2010 23:50:30
% Release Version 1.0
%
% ==========================================
if (nargin<6)
    mm = 0;
end
% Creates the MVAAR model
XDB = matrix_former(ts,p,mm);

% Total range over which to perform DTF
tot_range = [low_freq:high_freq];
% Number of frequencies
nfre = length(tot_range);

nchan = size(XDB, 1);
total = size(XDB, 4);
gamma2_set = zeros(total,nchan, nchan, nfre);
for i = 1:total
    TemA = cat(3, -eye(nchan), XDB(:, :, :, i));
    gamma2_set(i,:,:,:) = DTFvalue(TemA,low_freq,high_freq,fs);
end

end
