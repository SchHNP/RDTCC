function x = mvaar(y,p,UC)
% Multivariate (Vector) adaptive AR estimation base on a multidimensional
% Kalman filer algorithm. A standard VAR model (A0=I) is implemented. The
% state vector is defined as X=(A1|A2...|Ap) and x=vec(X')
%
% [x,e,Kalman,Q2] = mvaar(y,p,UC,mode,Kalman)
%
% The standard MVAR model is defined as:
%
%		y(n)-A1(n)*y(n-1)-...-Ap(n)*y(n-p)=e(n)
%
%	The dimension of y(n) equals s
%
%	Input Parameters:
%
% 		y			Observed data or signal
% 		p			prescribed maximum model order (default 1)
%		UC			update coefficient	(default 0.001)
%		mode	 	update method of the process noise covariance matrix 0...4 ^
%					correspond to S0...S4 (default 0)
%
%	Output Parameters
%
%		e			prediction error of dimension s
%		x			state vector of dimension s*s*p
%		Q2			measurement noise covariance matrix of dimension s x s
%		            测量噪声的协方差矩阵
%		
%

% Copyright (C) 2001-2002 Christian Kasess
%       $Revision: 1.3 $
%       $Id: mvaar.m,v 1.3 2005/05/25 13:02:03 schloegl Exp $
% Modifications (C) 2003 Alois Schloegl <a.schloegl@ieee.org>
%	docu improved
%	check for isnan(ERR) included
%	code straightened

%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% version 2 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.

if nargin<3,     UC = 0.001;      end
if nargin<2,     p = 1;           end
[M,LEN] = size(y');		%number of channels, total signal length
L = M*M*p;
if nargout>1,    x = zeros(L,LEN);        end

Kalman_Q2 = eye(M);
Kalman_Q1 = eye(L)*UC;
Kalman_Kp = eye(L);
Kalman_x = zeros(L,1);
% h = waitbar(0,'Begin to creat the mvar model, please wait...');
% pause(0.4)

for n = 2:LEN
    if(n<=p)
        Yr = [y(n-1:-1:1,:)' zeros(M,p-n+1)];	%vector of past observations
        Yr = Yr(:)';
    else
        Yr = y(n-1:-1:n-p,:)';						%vector of past observations
        Yr = Yr(:)';
    end
    %Update of measurement matrix
    Kalman_H = kron(eye(M),Yr);    
    %calculate prediction error
    ye = (Kalman_H*Kalman_x)';
    err = y(n,:)-ye;    
    if ~any(isnan(err(:)))
        %update of Q2 using the prediction error of the previous step
        Kalman_Q2 = (1-UC)*Kalman_Q2+UC*err'*err;        
        KpH = Kalman_Kp*Kalman_H';
        HKp = Kalman_H*Kalman_Kp;
        %Kalman gain
        Kalman_G = KpH/(Kalman_H*KpH+Kalman_Q2);
        %calculation of the a-posteriori state error covariance matrix
        %K=Kalman_Kp-Kalman_G*KpH'; Althouh PK is supposed to be symmetric, this operation makes the filter unstable
        K = Kalman_Kp-Kalman_G*HKp;
        %a-priori state error covariance matrix for the next time step
        Kalman_Kp = K+Kalman_Q1;
        %current estimation of state x
        Kalman_x = Kalman_x+Kalman_G*(err)';
    end % isnan>(err)
    x(:,n) = Kalman_x;
%     PerStr = fix((n/LEN)*100);
end

x = x';


