function [bic,aic] = mvar_model_order(X,MINP,MAXP)
%-----------------------------------------------------------------------
% FUNCTION: cca_find_model_order.m
% PURPOSE:  use the Bayesian Information Criterion (BIC) and/or the Aikaike
%           Information Criterion to find the best model order (NLAGS) 
%           for a multivariate data set.
% 
% INPUTS:   X: matrix of nvar variables by nobs observations of each variable    
%           MINP: minimum model order to consider
%           MAXP: maximum model order to consider
%
% OUTPUT:   bic: optimal model order according to BIC
%           aic: optimal model order according to Akaike Information Criterion (AIC)
%
%           Written by Anil Seth, March 2004
%           Updated August 2004
%           Updated December 2005
%           Ref: Seth, A.K. (2005) Network: Comp. Neural. Sys. 16(1):35-55
%-----------------------------------------------------------------------

[nvar,nobs] = size(X);
if(nobs<nvar);  error('Fewer observations than variables, exiting');    end
if(MAXP<=MINP); error('MAXP must be bigger than MINP, exiting');        end

bc = ones(1,MAXP).*999;
ac = ones(1,MAXP).*999;
for i = MINP:MAXP
%     eval('res = cca_regress_v1(X,i);','res = -1');   % estimate regression model, catch errors
    try
        res = cca_regress_v1(X,i);
    catch
        res = -1;
    end
    if(~isnumeric(res))
        [bc(i),ac(i)] = findbic(res,nvar,nobs,i);
        %disp(['VAR order ',num2str(i),', BIC = ',num2str(bc(i)),', AIC = ',num2str(ac(i))]);
    else
        disp('VAR failed');
        bc(i) = 999; 
        ac(i) = 999;
    end
end

[bicmin,bic] = min(bc);
[aicmin,aic] = min(ac);
end

%---------------------------------------------------------------------
function [bc,ac] = findbic(res,nvar,nobs,nlag)

error = log(abs(det(res.Z)));
nest = nvar*nvar*nlag;       
bc = error + (log(nobs)*nest/nobs);   
ac = error + (2*nest/nobs);
end

% -----------------------------------------------------------------------
function [ret] = cca_regress_v1(X, nlags)
% 
%   FUNCTION: cca_regress.m
%   PURPOSE:  perform multivariate regression
%
%   INPUT:  X           -   nvar (rows) by nobs (cols) observation matrix
%           NLAGS       -   number of lags to include in model
%
%   OUTPUT: ret.Z       -   covariance matrix of residuals

% figure regression parameters
[nvar,nobs] = size(X);
if(nvar>nobs); error('nvar>nobs, check input matrix'); end
% remove sample means if present (no constant terms in this regression)
% m = mean(X');
m = mean(X, 2);
if(abs(sum(m)) > 0.0001)
    mall = repmat(m,1,nobs);
    X = X-mall;
end
% if(abs(sum(m)) > 0.0001)
%     mall = repmat(m',1,nobs);
%     X = X-mall;
% end
% construct lag matrices
lags = -999*ones(nvar,nobs-nlags,nlags);
for jj=1:nvar
    for ii=1:nlags
        lags(jj,:,nlags-ii+1) = X(jj,ii:nobs-nlags+ii-1);
    end
end
%  regression (no constant term)
regressors = permute(lags, [2, 3, 1]);
regressors = reshape(regressors, size(regressors, 1), []);
% regressors = zeros(nobs-nlags,nvar*nlags);
% for ii=1:nvar
%     s1 = (ii-1)*nlags+1;
%     regressors(:,s1:s1+nlags-1) = squeeze(lags(ii,:,:));
% end

xdep = X(:, nlags+1:end)';
if rank(regressors) < min(size(regressors,1),size(regressors,2))
    beta = pinv(regressors)*xdep;
else
    beta = regressors\xdep;
end
xpred = regressors*beta;  % keep hold of predicted values
u = xdep-xpred;
ret.Z = cov(u);

% for ii=1:nvar
%     xvec = X(ii,:)';
%     xdep = xvec(nlags+1:end);
%     if rank(regressors) < ...
%             min(size(regressors,1),size(regressors,2))
% %         temp = regressors'*regressors;
% %         [u1,s,v] = svd(temp);
% %         s = s + 0.005*eye(size(s));
% %         temp = u1*s*v';
% %         beta(:,ii) = inv(temp)*regressors'*xdep;*/
%        beta(:,ii)=pinv(regressors)*xdep;
%     else
%        beta(:,ii) = regressors\xdep;
%     end  
% %     beta(:,ii) = regressors\xdep;
%     xpred(:,ii) = regressors*beta(:,ii);  % keep hold of predicted values
%     u(:,ii) = xdep-xpred(:,ii);
% end
%   organize output structure
% ret.Z = cov(u);
end






