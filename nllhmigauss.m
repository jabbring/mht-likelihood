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


