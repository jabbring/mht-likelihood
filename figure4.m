% //////////////////////////////////////////////////////////////////////
% Abbring and Salimans (2021), Figure 4 (fka strikes/mhte722I.src)
% - Reads preferred estimates, calculates Weibull MPH estimates, and 
%   writes data for Figure 4. Aggregate Strike End Hazard Rates
%
%
% dependencies: tab1.mat igausscdf.m igausspdf.m weibullcdf.m 
%               weibullpdf.m
% output:   fig4.csv - data aggregate strike end hazard rates
%           fig4.tex - LaTeX macros MPH results cited in text
%           weibullmph.mat - MPH results for gradient check (numgrad.m) 
% //////////////////////////////////////////////////////////////////////

%% clear screen and workspace and set seed
clear
clc
format short

rng(230676) % set seed for random starting values MPH estimation

%% settings
dispplot = true; % set to true to have script plot results

%% read strike data  
rawdata=load('strkdur.asc');
x=rawdata(:,2);
y=rawdata(:,1)/7;

nrobs=size(y,1);

%% calculate kernel smoothed empirical hazard rate
t=1/7:1/7:max(y);
rawcts=sum(abs(t-y)<1e-3)';
rawsrv=nrobs-[0;cumsum(rawcts(1:end-1))];
rawhaz=7*rawcts./rawsrv;

emphaz=7*ksdensity(t',t','Weights',rawhaz,'Kernel','epanechnikov',...
    'Bandwidth',1,'Support',[0,max(y)+1],'BoundaryCorrection','reflection');

%% calculate hazard rate estimated MHT model Table 1 IV
mht=load('tab1.mat'); % this loads the structure 'est' with MHT estimates
mu=1;
var=mht.est.bm_var;
p=mht.est.unobs_p;
v=mht.est.unobs_v;
exb=exp(x*mht.est.beta);

igpdf=0;
igcdf=0;
for i=1:nrobs
    igpdf=igpdf+igausspdf(t',v*exb(i)/mu,(v*exb(i)).^2/var)*p'/nrobs;
    igcdf=igcdf+igausscdf(t',v*exb(i)/mu,(v*exb(i)).^2/var)*p'/nrobs;
end
ighaz=igpdf./(1-igcdf);

%% estimate hazard rate Weibull MPH model
disp('Estimating Weibull MPH model');
options=optimset('GradObj','on','LargeScale','off','HessUpdate','bfgs',...
    'Display','off','MaxIter',1000);

% initialize
llh=[];
iter=0;
maxllh=-Inf;
pd=false;
% maximize likelihood with random starting points to find global optimum
while (length(llh)<3 || sum(llh>max(llh)-(1e-3)*nrobs)<3 || ~pd) && iter<10

    % new iteration
    iter=iter+1;
    disp(['Maximizing log likelihood, random initialization ' int2str(iter)]);
    
    startvalues = randn(1,9);
    [estpar,nllh,eflag,output,grad,hessian]=fminunc(@(par)nllhmph(par,y,...
        false,x,4),startvalues,options);

    % check and save results
    if eflag>0
        llh=[llh -nllh];
        if -nllh>maxllh
            maxllh=-nllh;
            par=estpar;

            % check positive definiteness of hessian
            try
                pd=all(diag(chol(hessian))>0);
            catch
                pd=false;
            end
        end
    end
end

% examine results
if ~isfinite(maxllh)
    error('All optimization attempts failed, check the model specification.')
elseif iter==10
    disp('Solution might be a local optimum!')
    disp('You may want increase the number of random starting values tried.')
else
    disp('Optimization was successful!')
end
if ~pd
    disp('Warning: Hessian at optimum was not positive definite')
end
    
% retrieve parameters
delta=exp(par(1));
v=exp(par(2:5));
p=exp([par(6:8) 0]);
p=p/sum(p);
beta=par(9);
exb=exp(x*beta);
wpdf=0;
wcdf=0;
for i=1:nrobs
    wpdf=wpdf+weibullpdf(t',exb(i).*v,delta)*p'/nrobs;
    wcdf=wcdf+weibullcdf(t',exb(i).*v,delta)*p'/nrobs;
end
weibullhaz=wpdf./(1-wcdf);

save('weibullmph','par'); % for checking numerical gradient

%% plot all three hazard rates if dissplot=true
if dispplot
    fig4 = figure('PaperSize',[20.98 29.68],'Color',[1 1 1]);
    axes4 = axes('Parent',fig4,'LineWidth',1,'FontSize',14);
    box(axes4,'on');
    hold(axes4,'all');
    plot(t,emphaz,'LineWidth',2,'Color','green');
    plot(t,weibullhaz,'LineWidth',2,'Color','blue');
    plot(t,ighaz,'LineWidth',2,'Color','red');
    title('Figure 4. Aggregate Strike End Hazard Rates');
    xlabel({'Strike duration in weeks'},'LineWidth',2,'FontSize',14);
    ylabel({'Hazard rate per week'},'LineWidth',2,'FontSize',14);
end

%% Export data to csv file for TikZ
f1=fopen('fig4.csv','w');   % mhtehazard.txt
fprintf(f1,'t, mht, mph, data, dummy\n');
fprintf(f1,'%6.6f, %6.6f, %6.6f, %6.6f, 1.0\n',[t' ighaz weibullhaz emphaz]');
fclose(f1);

%% Export MPH likelihood comparison to tex file
f1=fopen('fig4.tex','w');   
fprintf(f1,'\\def\\mphllh{$%6.1f$}\n',maxllh);
fprintf(f1,'\\def\\diffllh{$%6.1f$}\n',mht.llh-maxllh);
fclose(f1);

