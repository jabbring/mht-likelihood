% //////////////////////////////////////////////////////////////////////
% Abbring and Salimans (2021)
% - Simulates data once and estimates MHT model
%
% dependencies: simmht.m mhtmle.m
% //////////////////////////////////////////////////////////////////////

%% clear screen and workspace and set seed
clear
clc
format short

rng(230670) % set seed for random starting values MPH estimation

%% clear screen and workspace
clear
clc
format short

%% parameters
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
shocktype='point'
nrshocks=1
lambda=0.5
nu=-1
%shocktype='gamma'
%nrshocks=1
%lambda=1
%nu=[2; 1] %[nu; rho]  
beta=randn(k,1)/k
x=randn(n,k)/4;
ymax=randraw('exp',0.1, n); %10000; % to prevent infinite values / infinite computation time for defects
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
y=simmht(n,x,1,var,unobstype,v,p,shocktype,lambda,nu,beta,ymax);
cens=(y==ymax);
mean(cens)

%% Estimation
[parameters,standard_errors,max_log_likelihood,optimization]=...
    mhtmle(y,cens,x,unobstype,shocktype,nrunobs,nrshocks)
% last 4 arguments are: unobs_type, shock_type, nrunobs, nrshocks
% unobs_type and shock_type can be either 'gamma' or 'point'