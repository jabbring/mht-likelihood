function [lh, grad]=lhmigauss(par,y,cens,x,nrunobs)
% ///////////////////////////////////////////////////////////////////////
% Calculate negative loglikelihood mixed inverse Gaussian model
% Calculate analytical derivatives
% //////////////////////////////////////////////////////////////////////

% get parameters
mu=par(1);
var=par(2);
v=[1; par(3:1+nrunobs)];
p=[1-sum(par(2+nrunobs:2*nrunobs)); par(2+nrunobs:2*nrunobs)];

% check x, beta
if isempty(x)
    k=0;
    funx=ones(size(y));
else
    k=size(x,2);
    beta=par(end+1-k:end);
    funx=exp(x*beta);
end

lh=zeros(size(y));
for i=1:nrunobs
    % noncensored observations
    if sum(~cens)>0
        [igpdf{i}, pdfgrad{i}]=igausspdf(y(~cens),v(i)*funx(~cens)/mu,(v(i)*funx(~cens)).^2/var);
        lh(~cens)=lh(~cens)+p(i)*igpdf{i};
    end

    % random right censored observations
    if sum(cens)>0
        [igcdf{i}, cdfgrad{i}]=igausscdf(y(cens),v(i)*funx(cens)/mu,(v(i)*funx(cens)).^2/var);
        lh(cens)=lh(cens)+p(i)*(1-igcdf{i});
    end
end


