% //////////////////////////////////////////////////////////////////////
% Abbring and Salimans (2021), Table 1 (fka laplace/test.m)
% - Maximum Likelihood Estimates for Kennan’s (1985) Strike Duration 
%   Data
%
% dependencies: strkdur.asc mhtmle.m
% output:   tab1.tex - LaTeX version of Table 1
%           tab1times.tex - LaTeX code with comp times Table 1
%           tab1.mat - contains structure `est` with estimates and 
%                       corresponding llh column IV for figure4.m
% //////////////////////////////////////////////////////////////////////

%% clear screen and workspace and set seed
clear
clc
format short

rng(230670) % set seed for random starting values MHT estimation

%% read strike data   
rawdata=load('strkdur.asc');
x=rawdata(:,2);
y=rawdata(:,1)/7;

%% estimation
estimates=nan(6,14);
stderrors=estimates;
loglik=nan(6,1);
comptime = [];
% Columns I-VI
for i = 1:6
    fprintf('Calculating Table 1 Column %1d\n',i)
    L = min(i,5); % nrunobs
    Q = max(i-5,0); % nrshocks
    tic;
    [est,ses,llh,opt]=mhtmle(y,false,x,'point','point',L,Q);
    comptime=[comptime;toc];
    % last 4 arguments are: unobs_type, shock_type, nrunobs (L), nrshocks 
    % (Q); unobs_type and shock_type can be either 'gamma' or 'point'
    % disp(opt)
    if i==4
        save('tab1','est','llh') % save estimates IV for Figure 4
    end
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

%% Export tex file with computation times
f1=fopen('tab1times.tex','w'); 
fprintf(f1,'Computation times (in seconds):');
fprintf(f1,'\\begin{tabular}{cccccc}');
fprintf(f1,'I&II&III&IV&V&VI\\\\');
fprintf(f1,'$%4.1f$&$%4.1f$&$%4.1f$&$%4.1f$&$%4.1f$&$%4.1f$',comptime);
fprintf(f1,'\\end{tabular}\n');
fclose(f1);

%% Export tex file with Table 1 (incl macros estimates cited in main text)
f1=fopen('tab1.tex','w'); 
fprintf(f1,'\\def\\vone{$%6.1f$}\n',estimates(4,4+1));
fprintf(f1,'\\def\\vtwo{$%6.1f$}\n',estimates(4,4+2));
fprintf(f1,'\\def\\vthree{$%6.1f$}\n',estimates(4,4+3));
fprintf(f1,'\\def\\vfour{$%6.1f$}\n',estimates(4,4+4));

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
