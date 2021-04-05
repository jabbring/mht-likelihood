% //////////////////////////////////////////////////////////////////////
% Abbring and Salimans (2021), Table 1 Column I-V (fka laplace/test.m)
% - Comparison with MLE based on Exact Mixed Inverse Gaussian Likelihood
%
% dependencies: strkdur.asc migaussmle.m
% output:   tab1mig.tex - LaTeX version of Table 1 inverse Gaussian
%           tab1migtimes.tex - LaTeX code with comp times inv Gaussian
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

%% Estimation
estimates=nan(6,14);
stderrors=estimates;
loglik=nan(6,1);
comptime=loglik;
% Columns I-VI
for l = 1:5
    fprintf('Calculating Table 1 (BM) Column %1d\n',l)
    tic
    [est,cov,llh,opt]=migaussmle(y,y*0,x,l);
    comptime(l)=toc;
    loglik(l)=llh;
 
    mu=est(1);
    var=est(2)/mu^2;
    v=[1; est(3:1+l)]/mu;
    p=[1-sum(est(2+l:2*l)); est(2+l:2*l)];
    beta=est(end);

    estimates(l,1)=var; % sigma^2
    estimates(l,4)=beta; % beta
    [estimates(l,5:4+l),srtidx]=sort(v); % v
    estimates(l,10:9+l)=p(srtidx); % pi
    
    dvar = [-2*var/mu 1/mu^2];
    stderrors(l,1)=sqrt(dvar*cov(1:2,1:2)*dvar');
    stderrors(l,4)=sqrt(cov(end,end)); 
    dv = [-v [zeros(1,l-1);eye(l-1)]]/mu; 
    stdv = sqrt(diag(dv*cov([1 3:1+l],[1 3:1+l])*dv'));
    stderrors(l,5:4+l)=stdv(srtidx);
    dp = [-ones(1,l-1);eye(l-1)];
    stdp = sqrt(diag(dp*cov(2+l:2*l,2+l:2*l)*dp'));
    stderrors(l,10:9+l)=stdp(srtidx);
end

%% Export tex file with computation times
f1=fopen('tab1migtimes.tex','w'); 
fprintf(f1,'Computation times (analytical IG; in seconds):');
fprintf(f1,'\\begin{tabular}{cccccc}');
fprintf(f1,'I&II&III&IV&V&VI\\\\');
fprintf(f1,'$%4.1f$&$%4.1f$&$%4.1f$&$%4.1f$&$%4.1f$&$%4.1f$',comptime);
fprintf(f1,'\\end{tabular}\n');
fclose(f1);

%% Export tex file with Table 1 (incl macros estimates cited in main text)
f1=fopen('tab1mig.tex','w'); 
fprintf(f1,'\\def\\vone{$%6.1f$}\n',estimates(4,4+1));
fprintf(f1,'\\def\\vtwo{$%6.1f$}\n',estimates(4,4+2));
fprintf(f1,'\\def\\vthree{$%6.1f$}\n',estimates(4,4+3));
fprintf(f1,'\\def\\vfour{$%6.1f$}\n',estimates(4,4+4));

fprintf(f1,'%s\n','\begin{table}');
fprintf(f1,'%s\n','\caption{Replicating Table \ref{table:strike} Using Inverse Gaussian Pdf}');
fprintf(f1,'%s\n','\vspace*{0.5em}');
fprintf(f1,'%s\n','\begin{center}');
fprintf(f1,'%s\n','\small{\begin{tabular}{ccccccc}');
fprintf(f1,'%s\n','\toprule');
fprintf(f1,'%s\n','& I & II & III & IV & V & VI\tabularnewline');
fprintf(f1,'%s\n','\midrule');
fprintf(f1,'%s\n','\midrule');
fprintf(f1,'%s\n','$\mu$ & $1$ & $1$ & $1$ & $1$ & $1$ &  \tabularnewline');
fprintf(f1,'%s\n','& $(0)$ & $(0)$ & $(0)$ & $(0)$ & $(0)$ &  \tabularnewline');
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
    s=sprintf('$\\pi_%d$ & $%6.0f$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$ & $%6.3f$\\tabularnewline',...
        [l;estimates(:,9+l)]);
        fprintf(f1,'%s\n',strrep(s,'NaN',''));
    s=sprintf('& $(%6.0f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$ & $(%6.3f)$\\tabularnewline',...
        stderrors(:,9+l));
        fprintf(f1,'%s\n',strrep(s,'(   NaN)',''));
    fprintf(f1,'%s\n','\midrule');
end
fprintf(f1,'%s\n','\midrule');
s=sprintf('$\\ell_N$ & $%6.1f$ & $%6.1f$ & $%6.1f$ & $%6.1f$ & $%6.1f$ & $%6.1f$\\tabularnewline\n',...
    loglik);
fprintf(f1,'%s\n',strrep(s,'NaN',''));
fprintf(f1,'%s\n','\bottomrule');
fprintf(f1,'%s\n','\end{tabular}}');
fprintf(f1,'%s\n','\end{center}');
fprintf(f1,'%s\n','{\footnotesize Note: The drift is normalized to $1$ per week. All specifications include a single covariate, \cites{jem85:kennan} deseasonalized and detrended log industrial production.  Asymptotic standard errors are in parentheses.}');
fprintf(f1,'%s\n','\end{table}');
fclose(f1);

