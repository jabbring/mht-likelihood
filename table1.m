% //////////////////////////////////////////////////////////////////////
% Abbring and Salimans (2021), Table 1 (fka laplace/test.m)
% - Maximum Likelihood Estimates for Kennan’s (1985) Strike Duration 
%   Data
% //////////////////////////////////////////////////////////////////////

%% clear screen and workspace
clear
clc
format short

%% read strike data   
rawdata=load('strkdur.asc');
x=rawdata(:,2);
y=rawdata(:,1)/7;

%% estimation
estimates=nan(6,15);
stderrors=estimates;
loglik=nan(6,1);
estimates(:,1)=ones(6,1); % mu
stderrors(:,1)=zeros(6,1);
rng(230670); % seed for random start values
% Columns I-VI
for i = 1:6
    fprintf('Calculating Table 1 Column %1d\n',i)
    L = min(i,5); % nrunobs
    Q = max(i-5,0); % nrshocks
    [est,ses,llh,opt]=mhtmle(y,false,x,'point','point',L,Q);
    % last 4 arguments are: unobs_type, shock_type, nrunobs (L), nrshocks 
    % (Q); unobs_type and shock_type can be either 'gamma' or 'point'
    disp(opt)
    loglik(i)=llh;
    estimates(i,2)=est.bm_var; % sigma^2
    if Q>0
        estimates(i,3)=est.shock_lambda;
        estimates(i,4)=est.shock_nu;
    end
    estimates(i,5)=est.beta; % beta
    [estimates(i,6:5+L),srtidx]=sort(est.unobs_v); % v
    estimates(i,11:10+L)=est.unobs_p(srtidx); % pi
    
    stderrors(i,2)=ses.bm_var; 
    if Q>0
        stderrors(i,3)=ses.shock_lambda;
        stderrors(i,4)=ses.shock_nu;
    end
    stderrors(i,5)=ses.beta; 
    stderrors(i,6:5+L)=ses.unobs_v(srtidx);
    stderrors(i,11:10+L)=ses.unobs_p(srtidx);
end

%% Export tex file with Table 1

f1=fopen('table1.tex','w'); 
fprintf(f1,'logtlow, hist, dummy\n');
fprintf(f1,'%6.6f, %6.6f, 1.0\n',[xibounds' [fhist';NaN]]');
fclose(f1);