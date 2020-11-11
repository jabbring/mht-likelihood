% //////////////////////////////////////////////////////////////////////
% Abbring and Salimans (2021), Table 1 Column I-V (fka laplace/test.m)
% - Comparison with MLE based on Exact Mixed Inverse Gaussian Likelihood
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
y=rawdata(:,1)/7;

%% Estimation
tic
[parameters,standard_errors,max_log_likelihood,optimization]=...
    mhtmle(y,false,x,'point','point',1,1)
toc
% last 4 arguments are: unobs_type, shock_type, nrunobs, nrshocks
% unobs_type and shock_type can be either 'gamma' or 'point'


%% Estimation
tic
[parMIG, stderrMIG]=migaussmle(y,y*0,x,5)
toc