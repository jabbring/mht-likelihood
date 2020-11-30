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
estimates=nan(6,14);
stderrors=estimates;
loglik=nan(6,1);
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
    estimates(i,1)=est.bm_var; % sigma^2
    if Q>0
        estimates(i,2)=est.shock_lambda;
        estimates(i,3)=est.shock_nu;
    end
    estimates(i,4)=est.beta; % beta
    [estimates(i,5:4+L),srtidx]=sort(est.unobs_v); % v
    estimates(i,10:9+L)=est.unobs_p(srtidx); % pi
    
    stderrors(i,1)=ses.bm_var; 
    if Q>0
        stderrors(i,2)=ses.shock_lambda;
        stderrors(i,3)=ses.shock_nu;
    end
    stderrors(i,4)=ses.beta; 
    stderrors(i,5:4+L)=ses.unobs_v(srtidx);
    stderrors(i,10:9+L)=ses.unobs_p(srtidx);
end

%% Export tex file with Table 1

f1=fopen('tab1.tex','w'); 
fprintf(f1,'%s\n','\begin{table}');
fprintf(f1,'%s\n','\caption{Maximum Likelihood Estimates for \cites{jem85:kennan} Strike Duration Data\label{table:strike}}');
fprintf(f1,'%s\n','\vspace*{0.5em}');
fprintf(f1,'%s\n','\begin{center}');
fprintf(f1,'%s\n','\small{\begin{tabular}{ccccccc}');
fprintf(f1,'%s\n','\toprule');
fprintf(f1,'%s\n','& I & II & III & IV & V & VI\tabularnewline');
fprintf(f1,'%s\n','\midrule');
fprintf(f1,'%s\n','\midrule');
fprintf(f1,'%s\n','$\mu$ & $1$ & $1$ & $1$ & $1$ & $1$ & $1$ \tabularnewline');
fprintf(f1,'%s\n','& $(0)$ & $(0)$ & $(0)$ & $(0)$ & $(0)$ & $(0)$ \tabularnewline');
fprintf(f1,'%s\n','\midrule');
s=sprintf('$\\sigma^{2}$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$\\tabularnewline',...
    estimates(:,1));
    fprintf(f1,'%s\n',strrep(s,'NaN',''));
s=sprintf('& $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$\\tabularnewline',...
    stderrors(:,1));
    fprintf(f1,'%s\n',strrep(s,'(   NaN)',''));
fprintf(f1,'%s\n','\midrule');
s=sprintf('$\\lambda$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$\\tabularnewline',...
    estimates(:,2));
    fprintf(f1,'%s\n',strrep(s,'NaN',''));
s=sprintf('& $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$\\tabularnewline',...
    stderrors(:,2));
    fprintf(f1,'%s\n',strrep(s,'(   NaN)',''));
fprintf(f1,'%s\n','\midrule');
s=sprintf('$\\nu$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$\\tabularnewline',...
    estimates(:,3));
    fprintf(f1,'%s\n',strrep(s,'NaN',''));
s=sprintf('& $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$\\tabularnewline',...
    stderrors(:,3));
    fprintf(f1,'%s\n',strrep(s,'(   NaN)',''));
fprintf(f1,'%s\n','\midrule');
s=sprintf('$\\beta$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$\\tabularnewline',...
    estimates(:,4));
    fprintf(f1,'%s\n',strrep(s,'NaN',''));
s=sprintf('& $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$\\tabularnewline',...
    stderrors(:,4));
    fprintf(f1,'%s\n',strrep(s,'(   NaN)',''));
fprintf(f1,'%s\n','\midrule');
for l=1:5
    s=sprintf('$v_%d$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$\\tabularnewline',...
        [l;estimates(:,4+l)]);
        fprintf(f1,'%s\n',strrep(s,'NaN',''));
    s=sprintf('& $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$\\tabularnewline',...
        stderrors(:,4+l));
        fprintf(f1,'%s\n',strrep(s,'(   NaN)',''));
    fprintf(f1,'%s\n','\midrule');
end
for l=1:5
    s=sprintf('$\\pi_%d$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$\\tabularnewline',...
        [l;estimates(:,9+l)]);
        fprintf(f1,'%s\n',strrep(s,'NaN',''));
    s=sprintf('& $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$\\tabularnewline',...
        stderrors(:,9+l));
        fprintf(f1,'%s\n',strrep(s,'(   NaN)',''));
    fprintf(f1,'%s\n','\midrule');
end
fprintf(f1,'%s\n','\midrule');
fprintf(f1,'$\\ell_N$ & $%6.1f$ & $%6.1f$ & $%6.1f$ & $%6.1f$ & $%6.1f$ & $%6.1f$\\tabularnewline\n',...
    loglik);
fprintf(f1,'%s\n','\bottomrule');
fprintf(f1,'%s\n','\end{tabular}}');
fprintf(f1,'%s\n','\end{center}');
fprintf(f1,'%s\n','{\footnotesize Note: The drift is normalized to $1$ per week. All specifications include a single covariate, \cites{jem85:kennan} deseasonalized and detrended log industrial production.  Asymptotic standard errors are in parentheses.}');
fprintf(f1,'%s\n','\end{table}');
fclose(f1);