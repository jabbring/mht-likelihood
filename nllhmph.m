function [nllh, gradnllh]=nllhmph(par,y,cens,x,nrunobs)
% ///////////////////////////////////////////////////////////////////////
% Calculate negative loglikelihood mixed proportional hazards model
% //////////////////////////////////////////////////////////////////////

% censoring not implemented
if sum(cens)>0
        disp('WARNING: censoring not implemented in MPH model'); 
end

% get parameters
if size(par,1)>size(par,2);
    par=par';   % par is now row vector
end
delta=exp(par(1));
v=exp(par(2:1+nrunobs));
p=exp([par(2+nrunobs:2*nrunobs) 0]);
p=p/sum(p);

% check x, beta
if isempty(x)
    k=0;
    xb=ones(size(y));
else
    k=size(x,2);
    beta=par(end+1-k:end);
    exb=exp(x*beta);
end

clh=zeros(size(y,1),nrunobs);
grad=zeros(size(y,1),2*nrunobs+1);
for i=1:nrunobs
    [clh(:,i), gradi]=weibullpdf(y,exb.*v(i),delta);
    grad(:,1)=grad(:,1)+p(i)*delta*gradi(:,2);
    grad(:,1+i)=p(i).*exb.*v(i).*gradi(:,1);
    grad(:,2*nrunobs+1)=grad(:,2*nrunobs+1)+p(i).*x.*exb.*v(i).*gradi(:,1);
end
lh=clh*p';
grad(:,2+nrunobs:2*nrunobs)=...
    clh*([diag(p(1:nrunobs-1));zeros(1,nrunobs-1)]-p'.*p(1:nrunobs-1));

nllh=-sum(log(lh));
gradnllh = -sum(grad./lh);

