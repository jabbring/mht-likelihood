% //////////////////////////////////////////////////////////////////////
% Abbring and Salimans (2021), Figure 3 (fka laplace/PDF_comparison.m)
% - Approximate Probability Density and Histogram of Simulated Values 
%    of ln T for a Specification With Shocks and Heterogeneity
%
% dependencies: simmht.m numinvlap.m
% output:   fig3invlap.csv  - data approximate probability density
%           fig3hist.csv - data histogram
% //////////////////////////////////////////////////////////////////////

%% clear screen and workspace and set seed
clear
clc
format short

rng(230670);

%% settings
dispplot = false; % set to true to have script plot results

n=1e3; % sample size
neval = 250;
nbins = 100;
k=1; % number of regressors
var=1;

unobstype='point';
unobs_form=unobstype;
nrunobs=2;
v=[1; 5];
p=[0.7; 0.3];

shocktype='gamma';
unobs_shock=shocktype;
nrshocks=1;
lambda=1;
nu=[2; 1]; %[nu; rho] 

beta=0;

%% transformed parameter vector
par=log(var);
if isequal(unobstype,'point')
    par=[par; 1+log(p(2:end)/p(1)); log(v)];
else
    par=[par; log(v); log(p)];
end
par=[par; log(lambda); log(abs(nu)); beta];

%% generate data
x=ones(n,1);
ymax=1e10; % to prevent infinite values / infinite computation time for defects
% also to test censoring, can be scalar or vector
if isequal(shocktype,'gamma')
    eshock=lambda*nu(2)/nu(1);
else
    eshock=lambda'*nu;
end
if 1-abs(eshock)<=0
    warning('the chosen parameterization has a defect')
end

y=simmht(n,x,1,var,unobstype,v,p,shocktype,lambda,nu,beta,ymax);
cens=(y==ymax);

%mean_duration=mean(y)
%max_duration=max(y)

%% inverse LT

% smoothing
[fsmooth,xi] = ksdensity(log(y),'kernel','epanechnikov','width',...
        5*sqrt(mean(log(y).^2)-mean(log(y))^2)/sqrt(n),'npoints',nbins);
probs=numinvlap(eval(['@' unobstype shocktype]),par,exp(xi'),false,...
                                    ones(length(xi),1),nrunobs,nrshocks);

% histogram
if dispplot
    fig3 = figure('PaperSize',[20.98 29.68],'Color',[1 1 1]);
    intervalsize = (xi(end)-xi(1))/(nbins-1);
    xibounds = xi-intervalsize/2:intervalsize:xi(end)+intervalsize/2;
    nr=histc(log(y),xibounds); 
    fhist=nr(1:end-1)'/length(y)/intervalsize;
    axes3 = axes('Parent',fig3,'LineWidth',1,'FontSize',14);
    box(axes3,'on');
    hold(axes3,'all');
    bar((xibounds(2:end)+xibounds(1:end-1))/2,fhist,'FaceColor',[1 1 1],...
                                   'EdgeColor',[0 1 0],'BarWidth',0.5);
    plot(xi,probs,'LineWidth',2,'Color','blue');
    plot(xi,fsmooth,'LineWidth',2,'Color','red','LineStyle','--')
    title('Figure 3. Approximate PDF and Histogram of Simulated Values of ln T');
    xlabel({'ln t'},'LineWidth',2,'FontSize',14);
    ylabel({'Density of ln T'},'LineWidth',2,...
        'FontSize',14);
end

%% Export data to csv files for TikZ

f1=fopen('fig3invlap.csv','w'); % mcinvlap.csv
fprintf(f1,'logt, smooth, invlt, dummy\n');
fprintf(f1,'%6.6f, %6.6f, %6.6f, 1.0\n',[xi' fsmooth' probs]');
fclose(f1);

f1=fopen('fig3hist.csv','w'); % mchist.csv
fprintf(f1,'logtlow, hist, dummy\n');
fprintf(f1,'%6.6f, %6.6f, 1.0\n',[xibounds' [fhist';NaN]]');
fclose(f1);
