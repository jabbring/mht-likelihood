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

%% clear screen and workspace
clear
clc
format short

%% read strike data and preferred MHT estimates   
rawdata=load('strkdur.asc');
x=rawdata(:,2);
y=rawdata(:,1)/7;

nrobs=size(y,1);

%% empirical hazard rate
t=1/7:1/7:max(y);
rawcts=sum(abs(t-y)<1e-3)';
rawsrv=nrobs-[0;cumsum(rawcts(1:end-1))];
rawhaz=7*rawcts./rawsrv;

emphaz=7*ksdensity(t',t','Weights',rawhaz,'Kernel','epanechnikov',...
    'Bandwidth',1,'Support',[0,max(y)+1],'BoundaryCorrection','reflection');

%% hazard rate MHT model Table 1 IV

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

%% hazard rate Weibull MPH model
disp('Estimating Weibull MPH model');
rng(230676)
options=optimset('GradObj','on','LargeScale','off','HessUpdate','bfgs','Display','off','MaxIter',1000);

llh=[];
iter=0;
% maximize likelihood with random starting points to find global optimum
while (length(llh)<3 || sum(llh>max(llh)-(1e-3)*nrobs)<3) && iter<10
 
    % new iteration
    iter=iter+1;
    disp(['Maximizing log likelihood, random initialization ' int2str(iter)]);
    
    startvalues = randn(1,9);
    [estpar,nllh,eflag,output,grad,hessian]=fminunc(@(par)nllhmph(par,y,...
        false,x,4),startvalues,options);
    llh=[llh nllh];
end

delta=exp(estpar(1));
v=exp(estpar(2:5));
p=exp([estpar(6:8) 0]);
p=p/sum(p);
beta=estpar(9);
exb=exp(x*beta);
wpdf=0;
wcdf=0;
for i=1:nrobs
    wpdf=wpdf+weibullpdf(t',exb(i).*v,delta)*p'/nrobs;
    wcdf=wcdf+weibullcdf(t',exb(i).*v,delta)*p'/nrobs;
end
weibullhaz=wpdf./(1-wcdf);

save('weibullmph','estpar'); % for checking numerical gradient

figure(1)
plot(t,[emphaz weibullhaz ighaz]);

%% Export data to csv file for TikZ

f1=fopen('fig4.csv','w');   % mhtehazard.txt
fprintf(f1,'t, mht, mph, data, dummy\n');
fprintf(f1,'%6.6f, %6.6f, %6.6f, %6.6f, 1.0\n',[t' ighaz weibullhaz emphaz]');
fclose(f1);

%% Export MPH likelihood comparison to tex file

f1=fopen('fig4.tex','w');   
fprintf(f1,'\\def\\mphllh{$%6.1f$}\n',-nllh);
fprintf(f1,'\\def\\diffllh{$%6.1f$}\n',mht.llh+nllh);
fclose(f1);

