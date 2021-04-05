function [lapprod,lapgrad,Glapprod,Glapgrad]=pointpoint(par,s,y,cens,x,nrunobs,nrshocks)

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
if length(par)~=2*nrunobs+2*nrshocks+k
    error('wrong number of parameters')
end

% retrieve parameters
var=exp(par(1));
p=exp([1; par(2:nrunobs)]);
p=p/sum(p);
v=exp(par(nrunobs+1:2*nrunobs));
lambda=exp(par(2*nrunobs+1:2*nrunobs+nrshocks));
nu=-exp(par(2*nrunobs+nrshocks+1:2*nrunobs+2*nrshocks));
beta=par(end-k+1:end);

% laplace exponent without shocks
lapexp=s+0.5*var*s.^2;

% incorporate shocks
for j=1:nrshocks
    %lapexp=lapexp+lambda(j)*(exp(s*nu(j))-1-s*(nu(j)/(1+nu(j)^2)));
    lapexp=lapexp+lambda(j)*(exp(s*nu(j))-1); % JHA: removed compensation shocks
end

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
for j=1:nrshocks
    lapgrad=lapgrad+lambda(j)*nu(j)*exp(s*nu(j));
end
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
    for j=1:nrshocks
        Glapprod{2*nrunobs+j}=exp(loglapprod+log(exp(s*nu(j))-1)+logymat);
        Glapprod{2*nrunobs+nrshocks+j}=lambda(j)*exp(loglapprod+logymat+logs+s*nu(j));
    end
    
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
    for j=1:nrshocks
        Glapprod{2*nrunobs+j}=Glapprod{2*nrunobs+j}*lambda(j);
        Glapprod{2*nrunobs+nrshocks+j}=nu(j)*Glapprod{2*nrunobs+nrshocks+j};
    end
    
    % calculate derivative lapgrad without shocks
    Glapgrad{1}=s;
    Glapgrad{end}=var;
    
    % incorporate shocks
    for j=1:nrshocks
        Glapgrad{2*nrunobs+j}=nu(j)*exp(s*nu(j));
        Glapgrad{2*nrunobs+nrshocks+j}=lambda(j)*((1+nu(j)*s).*exp(s*nu(j)));
        Glapgrad{end}=Glapgrad{end}+lambda(j)*nu(j)^2*exp(s*nu(j));
    end
    
    % incorporate parametric structure
    Glapgrad{1}=Glapgrad{1}*var;
    for j=1:nrshocks
        Glapgrad{2*nrunobs+j}=Glapgrad{2*nrunobs+j}*lambda(j);
        Glapgrad{2*nrunobs+nrshocks+j}=nu(j)*Glapgrad{2*nrunobs+nrshocks+j};
    end
    
end