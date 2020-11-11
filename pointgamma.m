function [lapprod,lapgrad,Glapprod,Glapgrad]=pointgamma(par,s,y,cens,x,nrunobs,nrshocks)

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
if length(par)~=nrunobs*2+3+k
    error('wrong number of parameters')
end

% retrieve parameters
var=exp(par(1));
p=exp([1; par(2:nrunobs)]);
p=p/sum(p);
v=exp(par(nrunobs+1:2*nrunobs));
lambda=exp(par(2*nrunobs+1));
nu=exp(par(2*nrunobs+2));
rho=exp(par(2*nrunobs+3));
beta=par(end-k+1:end);

% laplace exponent
lapexp=s+0.5*var*s.^2;

% incorporate shocks
% lapexp=lapexp+lambda*(1./(s/nu+1).^rho-1+s*(rho*gammainc(nu,rho+1)/nu));
lapexp=lapexp+lambda*(1./(s/nu+1).^rho-1); % JHA: removed compensation shocks

% pre-allocate
lapprod=zeros(reqdim);

% xfun
if k>0
    logxfun=x*beta;
else
    logxfun=zeros(length(y),1);
end
xfun=exp(logxfun);

% calculate laplace transform at psi(s)
laparg=-s.*repmat(xfun,1,reqdim(2));
loglaparg=log(laparg);
logyprod=lapexp.*ymat+logymat.*repmat(~cens,1,reqdim(2));
for j=1:nrunobs
    lapprod=lapprod+p(j)*exp(laparg*v(j));
end

% multiply with exp(psi(s)*t)
loglapprod=log(lapprod)+logyprod;
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
    pgrad=cell(nrunobs,1);
    xfgrad=zeros(reqdim);
    
    % gradient laplace transform
    for j=1:nrunobs
        lad=laparg*v(j);
        if k>0
            xfgrad=xfgrad+v(j)*p(j)*exp(lad+loglaparg+logyprod);
        end
        if nrunobs>1
            pgrad{j}=exp(lad+logyprod);
        end
        Glapprod{nrunobs+j}=p(j)*exp(loglaparg+lad+logyprod);
        if ~isempty(Glapprod{end})
            Glapprod{end}=Glapprod{end}-p(j)*v(j)*exp(lad+repmat(logxfun,1,reqdim(2))+logyprod);
        else
            Glapprod{end}=-p(j)*v(j)*exp(lad+repmat(logxfun,1,reqdim(2))+logyprod);
        end
    end
    
    % correct for censoring
    if ~isempty(cnrs)
        xfgrad(cnrs,:)=-xfgrad(cnrs,:)./clapexp;
        for j=1:nrunobs
            if nrunobs>1
                pgrad{j}(cnrs,:)=-pgrad{j}(cnrs,:)./clapexp;
            end
            Glapprod{nrunobs+j}(cnrs,:)=-Glapprod{nrunobs+j}(cnrs,:)./clapexp;
        end
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
    Glapprod{2*nrunobs+1}=exp(loglapprod+log(1./(s/nu+1).^rho-1)+logymat);
    Glapprod{2*nrunobs+2}=(rho*s.*((s/nu+1).^(-rho-1))/nu^2).*exp(loglapprod+log(lambda)+logymat);
    Glapprod{2*nrunobs+3}=-log(s/nu+1).*exp(loglapprod+log(lambda)-log(s/nu+1)*rho+logymat);
    
    % incorporate parametric structure
    Glapprod{1}=Glapprod{1}*var;
    for j=1:k
        Glapprod{end-k-1+j}=repmat(x(:,j),1,reqdim(2)).*xfgrad;
    end
    for j=1:nrunobs
        Glapprod{nrunobs+j}=Glapprod{nrunobs+j}*v(j);
        ad=p(j)*pgrad{j};
        if j>1
            Glapprod{j}=Glapprod{j}+ad;
        end
        for l=2:nrunobs
            if ~isempty(Glapprod{l})
                Glapprod{l}=Glapprod{l}-p(l)*ad;
            else
                Glapprod{l}=-p(l)*ad;
            end
        end
    end
    Glapprod{2*nrunobs+1}=Glapprod{2*nrunobs+1}*lambda;
    Glapprod{2*nrunobs+2}=Glapprod{2*nrunobs+2}*nu;
    Glapprod{2*nrunobs+3}=Glapprod{2*nrunobs+3}*rho;
    
    % calculate derivative lapgrad without shocks
    Glapgrad{1}=s;
    Glapgrad{end}=var;
    
    % incorporate shocks
    Glapgrad{2*nrunobs+1}=rho*(-(s/nu+1).^(-rho-1))/nu;
    Glapgrad{2*nrunobs+2}=lambda*rho*((-rho-1)*s.*(s/nu+1).^(-rho-2)/nu^2+((s/nu+1).^(-rho-1))/nu)/nu;
    Glapgrad{2*nrunobs+3}=lambda*(-(s/nu+1).^(-rho-1)/nu+log(s/nu+1).*rho.*exp(log(s/nu+1)*(-rho-1))/nu);
    Glapgrad{end}=Glapgrad{end}+lambda*rho*(rho+1)*(s/nu+1).^(-rho-2)/nu^2;
    
    % incorporate parametric structure
    Glapgrad{1}=Glapgrad{1}*var;
    Glapgrad{2*nrunobs+1}=Glapgrad{2*nrunobs+1}*lambda;
    Glapgrad{2*nrunobs+2}=Glapgrad{2*nrunobs+2}*nu;
    Glapgrad{2*nrunobs+3}=Glapgrad{2*nrunobs+3}*rho;
    
end
