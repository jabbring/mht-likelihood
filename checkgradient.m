% //////////////////////////////////////////////////////////////////////
% Abbring and Salimans (2021), extra calculations (fka laplace/test3.m)
% - Check analytical gradients against numerical ones
%
% Dependencies: strkdur.asc weibullmph.mat numgrad.m mhtobj.m nllhmph.m
%               $(specs)
% Output: chckgrad.tex
% //////////////////////////////////////////////////////////////////////

%% clear screen and workspace
clear
clc
format short

rng(230670);

%% read strike data 
rawdata=load('strkdur.asc');
x=rawdata(:,2);
y=rawdata(:,1)/7;

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

[obj, pargrad] = mhtobj(eval(['@' unobs_form shock_form]),startvalues,y,...
    cens,x,nrunobs,nrshocks);
numgradient = numgrad(@(par)mhtobj(eval(['@' unobs_form shock_form]),...
    par,y,cens,x,nrunobs,nrshocks),startvalues);

f1=fopen('chckgrad.tex','w');   
fprintf(f1,'\\begin{table}[ht]\n');
fprintf(f1,'\\caption{Analytical and Numerical Gradients MHT}\n');
fprintf(f1,'\\begin{center}\n');
fprintf(f1,'\\begin{tabular}{rr}\n');
fprintf(f1,'$%6.6f$&$%6.6f$\\\\\n',[pargrad numgradient]');
fprintf(f1,'\\end{tabular}\n');
fprintf(f1,'\\end{center}\n');
fprintf(f1,'\\end{table}\n');

mph=load('weibullmph');

[obj,pargrad1] = nllhmph(mph.par,y,false,x,4);
numgradient1 = numgrad(@(par)nllhmph(par',y,false,x,4),mph.par');

altpar=mph.par.*(1+randn(size(mph.par))/10);
[obj,pargrad2] = nllhmph(altpar,y,false,x,4);
numgradient2 = numgrad(@(par)nllhmph(par',y,false,x,4),altpar');

fprintf(f1,'\\begin{table}[ht]\n');
fprintf(f1,'\\caption{Analytical and Numerical Gradients MPH}\n');
fprintf(f1,'\\begin{center}\n');
fprintf(f1,'\\begin{tabular}{rr|rr}\n');
fprintf(f1,'$%6.6f$&$%6.6f$&$%6.6f$&$%6.6f$\\\\\n',[pargrad1' numgradient1 pargrad2' numgradient2]');
fprintf(f1,'\\end{tabular}\n');
fprintf(f1,'\\end{center}\n');
fprintf(f1,'\\end{table}\n');
fclose(f1);