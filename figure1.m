% //////////////////////////////////////////////////////////////////////
% Abbring and Salimans (2021), Figure 1 (fka laplace/invtest.m)
% - Approximation Error of the Log Likelihood for Various M
%
% Dependencies: strkdur.asc migaussmle.m lhmigauss.m numinvlap.m 
%               numinvlap2.m pointpoint.m
% Output: fig1.csv fig1times.tex
% //////////////////////////////////////////////////////////////////////

%% clear screen and workspace and set seed
clear
clc
format long

est=false; % false=as in paper; true=use MIG ML estimates
rng(230670) % seed for random starting values

%% settings
dispplot = false; % set to true to have script plot results

nrunobs=4;
nrshocks=0;

nrsim = 100;
if est
        nrsim=1;
end
%% read strike data
rawdata=load('strkdur.asc');

x=rawdata(:,2);
y=rawdata(:,1);
cens=false(length(y),1);

timea = [];
timen = [];
for j=1:nrsim
fprintf('.');
%% parameters analytical, estimated or start

if est
    [parMLE,cov,llh,opt]=migaussmle(y,cens,x,nrunobs);
else
    % set starting values
    mm=mean(y);
    mu=1/mm;
    var=mean(1./y-mu);
    v=exp(randn(nrunobs-1,1));
    p=ones(nrunobs-1,1)/nrunobs;
    if isempty(x);
        beta=[];
    else
        beta=zeros(size(x,2),1);
    end
    parMLE=[mu; var; v; p; beta];
end

mu=parMLE(1);
var=parMLE(2);
v=[1; parMLE(3:1+nrunobs)];
p=[1-sum(parMLE(2+nrunobs:2*nrunobs)); parMLE(2+nrunobs:2*nrunobs)];
beta=parMLE(end);

%% parameters numerical

% rescale
var=var/mu^2;
v=v/mu;
p=log(p);
p=p-p(1)+1;
p=p(2:end);

% convert to parameter vector of startvalues
par2=[log(var); p; log(v); beta];

%% test
tic
a_probs=lhmigauss(parMLE,y,cens,x,nrunobs);
timea=[timea;toc];
tic
n_probs=numinvlap(@pointpoint,par2,y,cens,x,nrunobs,nrshocks)./y;
timen=[timen;toc];
errs=[a_probs abs(a_probs-n_probs) abs(a_probs-n_probs)./n_probs];
maxerrs=max(errs);
meanerrs=mean(abs(errs));
logerr=abs(sum(log(n_probs))-sum(log(a_probs)));

%% graph maker
le{j}=zeros(30,1);
for i=1:30
    le{j}(i)=log(abs(sum(log(numinvlap2(@pointpoint,par2,y,i,cens,x,...
        nrunobs,nrshocks)./y))-sum(log(a_probs))));
end
end
fprintf('\n');
lerr=zeros(30,1);
for j=1:nrsim
    lerr=lerr+le{j}/nrsim;
end
if dispplot
    fig1 = figure('PaperSize',[20.98 29.68],'Color',[1 1 1]);
    axes1 = axes('Parent',fig1,'YScale','log','YMinorTick','on','LineWidth',1,'FontSize',14);
    box(axes1,'on');
    hold(axes1,'all');
    plot(exp(lerr),'LineWidth',2,'Color','blue','Marker','o');
    title('Figure 1. Approximation Error of the Log Likelihood for Various M');
    xlabel({'M'},'LineWidth',2,'FontSize',14);
    ylabel({'Average absolute error in log likelihood'},'LineWidth',2,...
        'FontSize',14);
end

f1=fopen('fig1times.tex','w'); 
fprintf(f1,'Note: Mean calculation times for Figure 1 are $%0.5e$ seconds ',...
    mean(timea));
fprintf(f1,'(analytical) and $%0.5e$ seconds (numerical inversion), so that ',...
    mean(timen));
fprintf(f1,'mean time numerical $=%6.2f\\times$ mean time analytical.\n',...
    mean(timen)/mean(timea));
fclose(f1);

%% Export data to csv file for TikZ
f1=fopen('fig1.csv','w');                % mhtellherr.csv
fprintf(f1,'M, llherr, dummy\n');
fprintf(f1,'%6.0f, %6.6f, 1.0\n',[(1:30)' lerr/log(10)]');
fclose(f1);
