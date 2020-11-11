% //////////////////////////////////////////////////////////////////////
% Abbring and Salimans (2021), Table 1 (fka laplace/test3.m)
% - Max. Likelihood Estimates for Kennan’s (1985) Strike Duration Data 
% //////////////////////////////////////////////////////////////////////

%% clear screen and workspace
clear
clc
format short

%% read strike data
%if ispc;
%        rawdata=load('..\..\strikes\strkdur.asc');
%else
%        rawdata=load('../../strikes/strkdur.asc');
%end    
rawdata=load('strkdur.asc');
x=rawdata(:,2);
y=rawdata(:,1);

%% test
cens=false;
unobs_form='point';
shock_form='point';
nrunobs=5;
nrshocks=1;

% check input
if size(y,1)<size(y,2);
    y=y';
end
if size(x,1)<size(x,2);
    x=x';
end
if size(cens,1)<size(cens,2);
    cens=cens';
end
n=length(y);
k=size(x,2);
if ~all(size(y)==[n 1])
    error('input y has the wrong size')
end
if ~all(size(x)==[n k])
    error('input x has the wrong size')
end
cens=logical(cens);
if length(cens)==1
    cens=(ones(n,1)*cens>0);
end
if ~all(size(cens)==[n 1])
    error('input cens has the wrong size')
end

% initialize
llh=[];
iter=0;
maxllh=-Inf;

% scale of data
mm=mean(y);

% set starting values
var=log(mm^2*mean(1./y-1/mm))+randn/10;
if isequal(unobs_form,'point')
    v=log(ones(nrunobs,1)*mm)+randn(nrunobs,1)/2;
    p=ones(nrunobs-1,1);
elseif isequal(unobs_form,'gamma')
    if nrunobs>1
        error('nrunobs should be 1 for gamma specification')
    end
    v=randn/2;
    p=randn/2-log(mm);
else
    error('wrong distribution for unobservables')
end
k=size(x,2);
if k==0
    beta=[];
else
    beta=zeros(size(x,2),1);
end
if isequal(shock_form,'point')
    nu=log(mm)-log(4:1:3+nrshocks)'+randn(nrshocks,1)/2;
    lambda=-log(mm)*ones(nrshocks,1);
elseif isequal(shock_form,'gamma')
    if nrshocks>1
        error('nrshocks should be 1 for gamma specification')
    end
    lambda=-log(mm);
    nu=[randn/2; randn/2];
else
    error('wrong distribution for shocks')
end

% convert to parameter vector of startvalues
startvalues=[var; p; v; lambda; nu; beta];

[obj, pargrad] = mhtobj(eval(['@' unobs_form shock_form]),startvalues,y,cens,x,nrunobs,nrshocks)
numgrad(@(par)mhtobj(eval(['@' unobs_form shock_form]),par,y,cens,x,nrunobs,nrshocks),startvalues)