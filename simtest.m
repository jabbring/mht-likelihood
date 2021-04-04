% //////////////////////////////////////////////////////////////////////
% Abbring and Salimans (2021)
% - Simulates data once and estimates MHT model
%
% dependencies: simmht.m mhtmle.m migaussmle.m
% //////////////////////////////////////////////////////////////////////

%% clear screen and workspace and set seed
clear
clc
format short

rng(230670) 


%% parameters
M = 25; % numerical integration

n=10000 % sample size
k=1; % number of regressors

var=0.5

unobstype='point'
nrunobs=1
v=2
p=1
%unobstype='gamma'
%nrshocks=1
%v=1 %omega
%p=1 %tau

nrshocks=1 % set to 0 for inverse Gaussian model
lambda=1

shocktype='point'
nu=-0.5

%shocktype='gamma'
%nu=[2; 1] %[nu; rho]  

beta=randn(k,1)/k
%beta=0;

x=randn(n,k)/4;

%ymax=randraw('exp',0.5, n); %10000; % to prevent infinite values / infinite computation time for defects
ymax=10000; 
% also to test censoring, can be scalar or vector
if isequal(shocktype,'gamma')
    eshock=nu(2)/nu(1);
else
    eshock=abs(lambda'*nu);
end
if 1-eshock<0
    warning('the chosen parameterization has a defect')
end

%% generate data
if nrshocks==0
        shocktype='point'
        lambda=0
        nu=0
end
y=simmht(n,x,1,var,unobstype,v,p,shocktype,lambda,nu,beta,ymax);
cens=(y==ymax);
mean(cens)

%% Estimation
disp('Estimation using Laplace inversion');

[parameters,standard_errors,max_log_likelihood,optimization]=...
    mhtmle2(y,M,cens,x,unobstype,shocktype,nrunobs,nrshocks)
% last 4 arguments are: unobs_type, shock_type, nrunobs, nrshocks
% unobs_type and shock_type can be either 'gamma' or 'point'

%% Estimate using analytic form for inverse Gaussian likelihood
disp('Estimation Gaussian model using inverse Gaussian likelihood');

[par,cov,llh,opt]=migaussmle(y,cens,x,nrunobs)

mu=par(1);
par_mig=struct;
par_mig.bm_var=par(2)/mu^2;
par_mig.unobs_p=[1-sum(par(2+nrunobs:2*nrunobs)); par(2+nrunobs:2*nrunobs)];
par_mig.unobs_v=[1; par(3:1+nrunobs)]/mu;
par_mig.beta=par(end);

stde_mig=struct;
    dvar = [-2*par_mig.bm_var/mu 1/mu^2];
stde_mig.bm_var=sqrt(dvar*cov(1:2,1:2)*dvar');
    dp = [-ones(1,nrunobs-1);eye(nrunobs-1)];
stde_mig.unobs_p = sqrt(diag(dp*cov(2+nrunobs:2*nrunobs,2+nrunobs:2*nrunobs)*dp'));
    dv = [-par_mig.unobs_v [zeros(1,nrunobs-1);eye(nrunobs-1)]]/mu; 
stde_mig.unobs_v = sqrt(diag(dv*cov([1 3:1+nrunobs],[1 3:1+nrunobs])*dv'));
stde_mig.beta=sqrt(cov(end,end)); 

%% Output
disp('Estimates and standard errors')
parameters
par_mig
standard_errors
stde_mig


