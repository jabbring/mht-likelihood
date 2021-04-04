function [parMLE,cov,llh,opt]=migaussmle(y,cens,x,nrunobs)

% choose alternative information matrix estimator 
altim = 'fd'; % 'fd': Hessian calculated using finite diff analyt score
              % default: Hessian outputted by fminunc (BFGS)   

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
%options=optimset('Display','off','GradObj','on','AlwaysHonorConstraints','bounds','TolFun',tc,'TolX',tp); %'GradConstr','on',
options=optimset('Display','off','GradObj','on','AlwaysHonorConstraints','bounds');
[parMLE,nllh,opt]=fmincon(objfun,startvalues,A,b,[],[],lb,ub,[],options);

% obtain hessian at solution
if isequal(altim,'fd')
    disp('Using finite differences of analytical score to estimate information matrix!')
    hsc=1e-6;
    k=length(parMLE);
    hess=zeros(k,k);
    for i=1:k
        par1=parMLE;
        par1(i)=par1(i)-hsc;
        par2=parMLE;
        par2(i)=par2(i)+hsc;
        [nllh, grad1]=nllhmigauss(par1,y,cens,x,nrunobs);
        [nllh, grad2]=nllhmigauss(par2,y,cens,x,nrunobs);
        hess(:,i)=0.5*(grad2-grad1)/hsc;
    end
    hess=0.5*hess+0.5*hess';
else
  [parMLE,nllh,opt,output,grad,hess] = fminunc(objfun,...
                              parMLE+randn(size(parMLE))/100,options);
end
    
llh=-nllh;

% asymptotic standard errors
cov=inv(hess);
%stderr=sqrt(diag(cov));