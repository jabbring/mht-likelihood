function [parMLE, stderr]=migaussmle(y,cens,x,nrunobs)

% precision
tc=1e-6; % tolerance on negative loglikelihood value
tp=1e-6; % tolerance on parameter change

% scale of data
mm=mean(y);

% set starting values
mu=1/mm;
var=mean(1./y-mu);
v=exp(randn(nrunobs-1,1));
p=ones(nrunobs-1,1)/nrunobs;
if isempty(x);
    beta=[];
else
    beta=zeros(size(x,2),1);
end

% convert to parameter vector of startvalues
startvalues=[mu; var; v; p; beta];

% parameter bounds
lb=-Inf*ones(size(startvalues));
ub=Inf*ones(size(startvalues));
lb(2)=0; % variance >= 0
lb(3:1+nrunobs)=0; % v >= 0
lb(2+nrunobs:2*nrunobs)=0; % p >= 0

% inequality constraints (A*par<=b)
if ~isempty(p)
    A=[zeros(1,1+nrunobs) ones(1,nrunobs-1) zeros(1,size(x,2))]; % sum(p) <= 1
    b=1;
else
    A=[];
    b=[];
end

% maximize the loglikelihood function
objfun=@(par)nllhmigauss(par,y,cens,x,nrunobs);
ktropts=optimset('Display','iter','GradObj','on','AlwaysHonorConstraints','bounds','TolFun',tc,'TolX',tp); %'GradConstr','on',
%[parMLE,fval,flag]=ktrlink(objfun,startvalues,A,b,[],[],lb,ub,[],ktropts);
[parMLE,fval,flag]=fmincon(objfun,startvalues,A,b,[],[],lb,ub,[],ktropts);

% obtain hessian at solution
hsc=1e-6;
k=length(parMLE);
hess=zeros(k,k);
for i=1:k
    par1=parMLE;
    par1(i)=par1(i)-hsc;
    par2=parMLE;
    par2(i)=par2(i)+hsc;
    [llh, grad1]=nllhmigauss(par1,y,cens,x,nrunobs);
    [llh, grad2]=nllhmigauss(par2,y,cens,x,nrunobs);
    hess(:,i)=0.5*(grad2-grad1)/hsc;
end
hess=0.5*hess+0.5*hess';

% asymptotic standard errors
stderr=sqrt(diag(inv(hess)));