function [lapprod,lapgrad,Glapprod,Glapgrad]=gammagamma(par,s,y,cens,x,nrunobs,nrshocks)

% regressors?
k=size(x,2);

% resize censoring vector
if size(cens,1)<size(cens,2)
    cens=cens';
end
if length(cens)==1
    cens=logical(cens*ones(length(y),1));
end

% numbers censoring
hlpnrs=1:length(y);
cnrs=hlpnrs(cens);

% requested dimension
reqdim=size(s);
ymat=repmat(y,1,reqdim(2));
logymat=repmat(log(y),1,reqdim(2));
logs=log(s);

% check length of parameters
if length(par)~=6+k
    error('wrong number of parameters')
end

% retrieve parameters
var=exp(par(1));
omega=exp(par(2));
tau=exp(par(3));
lambda=exp(par(4));
nu=exp(par(5));
rho=exp(par(6));
beta=par(end-k+1:end);

% laplace exponent
lapexp=s+0.5*var*s.^2;

% incorporate shocks
% lapexp=lapexp+lambda*(1./(s/nu+1).^rho-1+s*(rho*gammainc(nu,rho+1)/nu)); 
lapexp=lapexp+lambda*(1./(s/nu+1).^rho-1); % JHA: removed compensation shocks

% xfun
if k>0
    logxfun=x*beta;
else
    logxfun=zeros(length(y),1);
end
xfun=exp(logxfun);

% calculate laplace transform at psi(s)
laparg=s.*repmat(xfun,1,reqdim(2));
loglaparg=log(laparg);
logyprod=lapexp.*ymat+logymat.*repmat(~cens,1,reqdim(2));
loglapprod=-tau*log(laparg/omega+1);

% multiply with exp(psi(s)*t)
loglapprod=loglapprod+logyprod;
lapprod=exp(loglapprod);

% correct for censoring
if ~isempty(cnrs)
    clapexp=lapexp(cnrs,:);
    logclapexp=log(clapexp);
    lapprod(cnrs,:)=exp(logyprod(cnrs,:)-logclapexp)-exp(loglapprod(cnrs,:)-logclapexp);
    loglapprod(cnrs,:)=log(lapprod(cnrs,:));
end

% calculate derivative of laplace exponent to s without shocks
lapgrad=1+var*s;

% incorporate shocks
lapgrad=lapgrad-lambda*rho*(s/nu+1).^(-rho-1)/nu;
loglapgrad=log(lapgrad);

% derivatives
if nargout==4
    
    % pre-allocate derivatives
    Glapprod=cell(length(par)+1,1);
    Glapgrad=cell(length(par)+1,1);
    
    % gradient laplace transform
    if k>0
    	xfgrad=-exp(log(tau/omega)+logs+logyprod+(-tau-1)*log((s.*repmat(xfun,1,reqdim(2)))/omega+1));
    end
    Glapprod{2}=exp(log(tau)-2*log(omega)+loglaparg+(-tau-1)*log(laparg/omega+1)+logyprod);
    Glapprod{3}=-log(laparg/omega+1).*lapprod;
    Glapprod{end}=-exp(log(tau/omega)+(-tau-1)*log(laparg/omega+1)+logyprod+repmat(logxfun,1,reqdim(2)));
    
    % correct for censoring
    if ~isempty(cnrs)
        xfgrad(cnrs,:)=-xfgrad(cnrs,:)./clapexp;
        Glapprod{2}(cnrs,:)=-Glapprod{2}(cnrs,:)./clapexp;
        Glapprod{3}(cnrs,:)=-Glapprod{3}(cnrs,:)./clapexp;
        Glapprod{end}(cnrs,:)=-Glapprod{end}(cnrs,:)./clapexp;
    end
    
    % incorporate gradients of laplace exponent
    if ~isempty(cnrs)
        % adjust ymat as short-cut to correct for censoring
        ymat(cnrs,:)=ymat(cnrs,:)-1./clapexp;
        logymat(cnrs,:)=log(ymat(cnrs,:));
    end    
    Glapprod{1}=0.5*exp(loglapprod+2*logs+logymat);
    Glapprod{end}=Glapprod{end}+exp(loglapprod+loglapgrad+logymat);
    Glapprod{4}=(1./(s/nu+1).^rho-1).*exp(loglapprod+logymat);
    Glapprod{5}=(rho*s.*((s/nu+1).^(-rho-1))/nu^2).*exp(loglapprod+log(lambda)+logymat);
    Glapprod{6}=-log(s/nu+1).*exp(loglapprod+log(lambda)-log(s/nu+1)*rho+logymat);

    % incorporate parametric structure
    Glapprod{1}=Glapprod{1}*var;
    for j=1:k
        Glapprod{end-k-1+j}=repmat(x(:,j),1,reqdim(2)).*xfgrad;
    end
    Glapprod{2}=Glapprod{2}*omega;
    Glapprod{3}=Glapprod{3}*tau;
    Glapprod{4}=Glapprod{4}*lambda;
    Glapprod{5}=Glapprod{5}*nu;
    Glapprod{6}=Glapprod{6}*rho;
    
    % calculate derivative lapgrad without shocks
    Glapgrad{1}=s;
    Glapgrad{end}=var;
    
    % incorporate shocks
    Glapgrad{4}=rho*(-(s/nu+1).^(-rho-1))/nu;
    Glapgrad{5}=lambda*rho*((-rho-1)*s.*(s/nu+1).^(-rho-2)/nu^2+((s/nu+1).^(-rho-1))/nu)/nu;
    Glapgrad{6}=lambda*(-(s/nu+1).^(-rho-1)/nu+log(s/nu+1).*rho.*exp(log(s/nu+1)*(-rho-1))/nu);
    Glapgrad{end}=Glapgrad{end}+lambda*rho*(rho+1)*(s/nu+1).^(-rho-2)/nu^2;
    
    % incorporate parametric structure
    Glapgrad{1}=Glapgrad{1}*var;
    Glapgrad{4}=Glapgrad{4}*lambda;
    Glapgrad{5}=Glapgrad{5}*nu;
    Glapgrad{6}=Glapgrad{6}*rho;
    
end
