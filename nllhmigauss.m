function [nllh, grad]=nllhmigauss(par,y,cens,x,nrunobs)
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

% negative loglikelihood
nllh=-sum(log(lh));

% gradients
if nargout==2
    grad=zeros(size(par));
    pgrad=zeros(nrunobs,1);
    funxgrad=zeros(size(y));
    for i=1:nrunobs
        % noncensored observations
        if sum(~cens)>0
            grad(1)=grad(1)+p(i)*sum(pdfgrad{i}(:,1).*(v(i)*funx(~cens)/mu^2)./lh(~cens));
            grad(2)=grad(2)+p(i)*sum(pdfgrad{i}(:,2).*((v(i)*funx(~cens)).^2/var^2)./lh(~cens));
            if i>1
                grad(i+1)=grad(i+1)-p(i)*sum((pdfgrad{i}(:,1).*funx(~cens)/mu+2*v(i)*pdfgrad{i}(:,2).*funx(~cens).^2/var)./lh(~cens));
            end
            pgrad(i)=pgrad(i)-sum(igpdf{i}./lh(~cens));
            funxgrad(~cens)=funxgrad(~cens)-p(i)*((pdfgrad{i}(:,1).*v(i)./mu+2*v(i)^2*pdfgrad{i}(:,2).*funx(~cens)/var)./lh(~cens));
        end

        % random right censored observations
        if sum(cens)>0
            grad(1)=grad(1)-p(i)*sum(cdfgrad{i}(:,1).*(v(i)*funx(cens)/mu^2)./lh(cens));
            grad(2)=grad(2)-p(i)*sum(cdfgrad{i}(:,2).*((v(i)*funx(cens)).^2/var^2)./lh(cens));
            if i>1
                grad(i+1)=grad(i+1)+p(i)*sum((cdfgrad{i}(:,1).*funx(cens)/mu+2*v(i)*cdfgrad{i}(:,2).*funx(cens).^2/var)./lh(cens));
            end
            pgrad(i)=pgrad(i)+sum(igcdf{i}./lh(cens));
            funxgrad(cens)=funxgrad(cens)+p(i)*((cdfgrad{i}(:,1).*v(i)./mu+2*v(i)^2*cdfgrad{i}(:,2).*funx(cens)/var)./lh(cens));
        end
    end
    if nrunobs>1
        grad(2+nrunobs:2*nrunobs)=pgrad(2:end)-pgrad(1);
    end
    for i=1:k
        grad(end-k+i)=(funx.*x(:,i))'*funxgrad;
    end
end

%function [pd, grad]=igausspdf(y,mu,lambda)
%pd=sqrt(lambda./(2*pi*(y.^3))).*exp(-lambda.*((y-mu).^2)./((2*mu.^2).*y));
%if nargout==2
%    grad=[pd.*(2*lambda.*(y-mu)./(2*mu.^2.*y)+4*mu.*y.*lambda.*(y-mu).^2./(2*mu.^2.*y).^2) ...
%        pd.*(0.5./lambda-(y-mu).^2./(2*mu.^2.*y))];
%end

%function [cd, grad]=igausscdf(y,mu,lambda)
%cd=0.5*erfc(-sqrt(lambda./y).*(y./mu-1)/sqrt(2))+exp(2*lambda./mu).*0.5.*erfc(sqrt(lambda./y).*(y./mu+1)./sqrt(2));
%if nargout==2
%    grad=[-exp(-(lambda./y).*(y./mu-1).^2/2)./sqrt(pi).*(sqrt(lambda./y).*(y./mu.^2)/sqrt(2))-...
%        exp(2*lambda./mu).*(lambda./mu.^2).*erfc(sqrt(lambda./y).*(y./mu+1)./sqrt(2))+...
%        exp(2*lambda./mu).*exp(-(lambda./y).*(y./mu+1).^2./2).*(sqrt(lambda./y).*(y./mu.^2)./sqrt(2))./sqrt(pi) ...
%        exp(-(lambda./y).*(y./mu-1).^2/2)./sqrt(pi).*(0.5*sqrt(y./lambda).*(1./y).*(y./mu-1)/sqrt(2))+...
%        exp(2*lambda./mu).*(1./mu).*erfc(sqrt(lambda./y).*(y./mu+1)./sqrt(2))-...
%        exp(2*lambda./mu).*exp(-(lambda./y).*(y./mu+1).^2./2).*(0.5*sqrt(y./lambda).*(1./y).*(y./mu+1)./sqrt(2))/sqrt(pi)];
%end



