function [par_out,stde_out,llh_out,opt_out] = mhtmle(y,cens,x,unobs_form,shock_form,nrunobs,nrshocks)
% ////////////////////////////////////////////////////////////////////////
% Maximum likelihood estimation of MHT model using numerical inversion of
% the laplace transform
%
% Author: Tim Salimans (minor adaptations by Jaap Abbring)
%
% INPUT:
%   - y, n x 1 double or integer vector of observations
%   - cens, n x 1 logical vector of scalar indicating whether the
%       observations are censored (1 for censoring, 0 not)
%   - x, n x k double matrix of explanatory variables
%   - unobs_form, string indicating the distribution of the unobservables
%       currently supports 'point' and 'gamma'
%   - shock_form, string indicating the distribution of the shocks
%       currently supports 'point' and 'gamma'
%   - nrunobs, scalar indicating the number of unobservables
%       only used for pointwise unobservables
%   - nrshocks, scalar indicating the number of shocks
%       only used for pointwise shocks
%
% OUTPUT:
%   - par_out, structure of estimated parameters
%   - stde_out, structure containing the asymptotic standard errors
%   - llh_out, scalar containing the maximum log likelihood value
%   - opt_out, structure containing information on the optimization
%
%   Information on convergence is displayed on screen.
%
% ////////////////////////////////////////////////////////////////////////

% choose alternative information matrix estimator 
altim = 'fd'; % 'fd': Hessian calculated using finite diff analyt score
              % 'op': OPG estimator
              % default: Hessian outputted by fminunc (BFGS)   
               
% check input
if size(y,1)<size(y,2);
    y=y';
end
if size(x,1)<size(x,2);
    x=x';
end
if size(cens,1)<size(cens,2);
    cens=cens';
end
n=length(y);
k=size(x,2);
if ~all(size(y)==[n 1])
    error('input y has the wrong size')
end
if ~all(size(x)==[n k])
    error('input x has the wrong size')
end
cens=logical(cens);
if length(cens)==1
    cens=(ones(n,1)*cens>0);
end
if ~all(size(cens)==[n 1])
    error('input cens has the wrong size')
end
if nargin<7
    if nargin<5
        error('too few input arguments')
    elseif nargin==5
        nrunobs=1;
        nrshocks=1;
    else
        if isequal(unobs_form,'gamma')&&isequal(shock_form,'point')
            nrshocks=nrunobs;
            nrunobs=1;
        elseif isequal(unobs_form,'point')&&isequal(shock_form,'gamma')
            nrshocks=1;
        else
            error('too few input arguments')
        end
    end
end

% initialize
llh=[];
iter=0;
maxllh=-Inf;
pd=false;

% maximize likelihood with random starting points to find global optimum
while (length(llh)<3 || sum(llh>max(llh)-(1e-3)*n)<3 || ~pd) && iter<10
    
    % new iteration
    iter=iter+1;
    disp(['Maximizing log likelihood, random initialization ' int2str(iter)]);
    
    % parameter initialization
    
    % scale of data
    mm=mean(y);
    
    % set starting values
    var=log(mm^2*mean(1./y-1/mm))+randn/10;
    if isequal(unobs_form,'point')
        v=log(ones(nrunobs,1)*mm)+randn(nrunobs,1)/2;
        p=ones(nrunobs-1,1);
    elseif isequal(unobs_form,'gamma')
        if nrunobs>1
            error('nrunobs should be 1 for gamma specification')
        end
        v=randn/2;
        p=randn/2-log(mm);
    else
        error('wrong distribution for unobservables')
    end
    k=size(x,2);
    if k==0
        beta=[];
    else
        beta=zeros(size(x,2),1);
    end
    if isequal(shock_form,'point')
        nu=log(mm)-log(4:1:3+nrshocks)'+randn(nrshocks,1)/2;
        lambda=-log(mm)*ones(nrshocks,1);
    elseif isequal(shock_form,'gamma')
        if nrshocks>1
            error('nrshocks should be 1 for gamma specification')
        end
        lambda=-log(mm);
        nu=[randn/2; randn/2];
    else
        error('wrong distribution for shocks')
    end
    
    % convert to parameter vector of startvalues
    startvalues=[var; p; v; lambda; nu; beta];
    
    % estimation
    options=optimset('GradObj','on','LargeScale','off','HessUpdate','bfgs','Display','off','MaxIter',1000);
    eval(['[estpar,nllh,eflag,output,grad,hessian]=fminunc(@(par)mhtobj(@'...
        unobs_form shock_form ',par,y,cens,x,nrunobs,nrshocks),startvalues,options);']);
            
    % check and save results
    if eflag>0
        llh=[llh -nllh];
        if -nllh>maxllh
            maxllh=-nllh;
            par=estpar;
            mlehess=hessian;
            optinf=output;
            
            % check positive definiteness of hessian
            try
                pd=all(diag(chol(mlehess))>0);
            catch
                pd=false;
            end
        end
    end
    
end

% examine results
if ~isfinite(maxllh)
    error('All optimization attempts failed, check the model specification.')
elseif iter==10
    disp('Solution might be a local optimum!')
    disp('Try using a simpler specification.')
else
    disp('Optimization was successful!')
end
if ~pd
    disp('Hessian at optimum was not positive definite, so standard errors are inaccurate.')
end

% retrieve parameters
par_out=struct;
par_out.bm_var=exp(par(1));
if isequal(unobs_form,'point')
    par_out.unobs_p=exp([1; par(2:nrunobs)])';
    par_out.unobs_p=par_out.unobs_p/sum(par_out.unobs_p);
    par_out.unobs_v=exp(par(nrunobs+1:2*nrunobs))';
    if nrshocks>0
        if isequal(shock_form,'point')
            par_out.shock_lambda=exp(par(2*nrunobs+1:2*nrunobs+nrshocks))';
            par_out.shock_nu=-exp(par(2*nrunobs+nrshocks+1:2*nrunobs+2*nrshocks))';
            shockind=(2*nrunobs+1:2*nrunobs+2*nrshocks);
        else
            par_out.shock_lambda=exp(par(2*nrunobs+1));
            par_out.shock_nu=exp(par(2*nrunobs+2));
            par_out.shock_rho=exp(par(2*nrunobs+3));
            shockind=(2*nrunobs+1:2*nrunobs+3);
        end
    end
else
    par_out.unobs_omega=exp(par(2));
    par_out.unobs_tau=exp(par(3));
    if nrshocks>0
        if isequal(shock_form,'point')
            par_out.shock_lambda=exp(par(4:3+nrshocks))';
            par_out.shock_nu=-exp(par(4+nrshocks:3+2*nrshocks))';
            shockind=(4:3+2*nrshocks);
        else
            par_out.shock_lambda=exp(par(4));
            par_out.shock_nu=exp(par(5));
            par_out.shock_rho=exp(par(6));
            shockind=(4:6);
        end
    end
end
par_out.beta=par(end-k+1:end)';

% alternative IM estimates
if isequal(altim,'fd')
    disp('Using finite differences of analytical score to estimate information matrix!')
    eval(['mlehess=numjac(@(dpar)mhtgrad(@'...
    unobs_form shock_form ',dpar,y,cens,x,nrunobs,nrshocks),par);']);
    mlehess=0.5*(mlehess+mlehess');
elseif isequal(altim,'op')
    disp('Using OPG estimator of information matrix!')
    eval(['[opprobs,opgrad]=numinvlap(@'...
        unobs_form shock_form ',par,y,cens,x,nrunobs,nrshocks);']);
    opgrad=-opgrad.*(1./opprobs);
    mlehess = opgrad'*opgrad;
end

% retieve asymptotic standard errors using delta rule
if nargout>=2
    
    % check whether all shock-parameters are identified
    if nrshocks>0
        if isequal(shock_form,'point')
            shockscale=-par_out.shock_lambda.*par_out.shock_nu;
        else
            shockscale=par_out.shock_lambda*par_out.shock_rho/par_out.shock_nu;
        end
        if all(shockscale>0.05)
            shocknotid=false;
            mlevar=inv(mlehess);
        else
            shocknotid=true;
            ind=1:length(par);
            ind(shockind)=[];
            mlevar=zeros(length(par),length(par));
            mlevar(ind,ind)=inv(mlehess(ind,ind));
        end
    else
        shocknotid=false;
        mlevar=inv(mlehess);
    end
        
    % get standard errors
    stde_out=struct;
    stde_out.bm_var=par_out.bm_var*sqrt(mlevar(1,1));
    if isequal(unobs_form,'point')
        gmat=[zeros(1,nrunobs-1); diag(par_out.unobs_p(2:end))]-par_out.unobs_p'*par_out.unobs_p(2:end);
        stde_out.unobs_p=sqrt(diag(gmat*mlevar(2:nrunobs,2:nrunobs)*gmat'))';
        stde_out.unobs_v=par_out.unobs_v.*sqrt(diag(mlevar(nrunobs+1:2*nrunobs,nrunobs+1:2*nrunobs)))';
        if nrshocks>0
            if isequal(shock_form,'point')
                if shocknotid
                    stde_out.shock_lambda='not identified';
                    stde_out.shock_nu='not identified';
                else
                    stde_out.shock_lambda=par_out.shock_lambda.*...
                        sqrt(diag(mlevar(2*nrunobs+1:2*nrunobs+nrshocks,2*nrunobs+1:2*nrunobs+nrshocks)))';
                    stde_out.shock_nu=-par_out.shock_nu.*...
                        sqrt(diag(mlevar(2*nrunobs+nrshocks+1:2*nrunobs+2*nrshocks,2*nrunobs+nrshocks+1:2*nrunobs+2*nrshocks)))';
                end
            else
                if shocknotid
                    stde_out.shock_lambda='not identified';
                    stde_out.shock_nu='not identified';
                    stde_out.shock_rho='not identified';
                else
                    stde_out.shock_lambda=par_out.shock_lambda*sqrt(mlevar(2*nrunobs+1,2*nrunobs+1));
                    stde_out.shock_nu=par_out.shock_nu*sqrt(mlevar(2*nrunobs+2,2*nrunobs+2));
                    stde_out.shock_rho=par_out.shock_rho*sqrt(mlevar(2*nrunobs+3,2*nrunobs+3));
                end
            end
        end
    else
        stde_out.unobs_omega=par_out.unobs_omega*sqrt(mlevar(2,2));
        stde_out.unobs_tau=par_out.unobs_tau*sqrt(mlevar(3,3));
        if nrshocks>0
            if isequal(shock_form,'point')
                if shocknotid
                    stde_out.shock_lambda='not identified';
                    stde_out.shock_nu='not identified';
                else
                    stde_out.shock_lambda=par_out.shock_lambda.*sqrt(diag(mlevar(4:3+nrshocks,4:3+nrshocks)))';
                    stde_out.shock_nu=-par_out.shock_nu.*sqrt(diag(mlevar(4+nrshocks:3+2*nrshocks,4+nrshocks:3+2*nrshocks)))';
                end
            else
                if shocknotid
                    stde_out.shock_lambda='not identified';
                    stde_out.shock_nu='not identified';
                    stde_out.shock_rho='not identified';
                else
                    stde_out.shock_lambda=par_out.shock_lambda*sqrt(mlevar(4,4));
                    stde_out.shock_nu=par_out.shock_nu*sqrt(mlevar(5,5));
                    stde_out.shock_rho=par_out.shock_rho*sqrt(mlevar(6,6));
                end
            end
        end
    end
    stde_out.beta=sqrt(diag(mlevar(end-k+1:end,end-k+1:end)))';
    
end

% final log likelihood
if nargout>=3
    llh_out=maxllh;
end

% optimization information
if nargout==4
    opt_out=optinf;
end
