% //////////////////////////////////////////////////////////////////////
% Abbring and Salimans (2021), Figure 2 (fka laplace/test2.m )
% - Approximation Error of the Log Inverse Gaussian Density Function
%
% Dependencies: numinvlap2 scatterplot
% Output: - fig2.csv
% //////////////////////////////////////////////////////////////////////

%% clear screen and workspace
clear
clc
format long

%% settings
dispplot = false; % set to true to have script plot results

%% calculate approximations
y=(0.1:0.1:44)';
lnp=log(numinvlap2(@pointpoint,[0; 0],y,25,false,[],1,0)./y);
lap=-0.5*log(2*pi)-1.5*log(y)-(y-1).^2./(2.*y);
err=abs(lnp-lap);
if dispplot
    scatterplot(lap,err)
end


%% Export data to csv file for TikZ

f1=fopen('fig2.csv','w');   % mhteproberr.csv
fprintf(f1,'lap, logerr, dummy\n');
fprintf(f1,'%6.6f, %6.6f, 1.0\n',[lap log10(err)]');
fclose(f1);