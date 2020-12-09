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
dispplot = true; % set to true to have script plot results

%% calculate approximations
y=(0.1:0.1:44)';
lnp=log(numinvlap2(@pointpoint,[0; 0],y,25,false,[],1,0)./y);
lap=-0.5*log(2*pi)-1.5*log(y)-(y-1).^2./(2.*y);
err=abs(lnp-lap);
if dispplot
    fig2 = figure('Color',[1 1 1]);
    axes2 = axes('Parent',fig2,'YScale','log','YMinorTick','on','LineWidth',1,...
        'FontSize',14);
    xlim(axes2,[-28 0]);
    ylim(axes2,[1e-12 1.5e1]);
    box(axes2,'on');
    hold(axes2,'all');
    scatter(lap,err,'MarkerFaceColor','blue','MarkerEdgeColor','blue',...
        'Marker','.');
    title('Figure 2. Approximation Error of the Log Inverse Gaussian PDF');
    xlabel('Log probability density','LineWidth',2,'FontSize',14);
    ylabel('Absolute error in log probability density','LineWidth',2,'FontSize',14);
end

%% Export data to csv file for TikZ

f1=fopen('fig2.csv','w');   % mhteproberr.csv
fprintf(f1,'lap, logerr, dummy\n');
fprintf(f1,'%6.6f, %6.6f, 1.0\n',[lap log10(err)]');
fclose(f1);