function grad=numgrad(functionhandle,par)
grad=zeros(size(par));
edif=1e-6;
for i=1:length(par)
    par1=par;
    par2=par;
    par1(i)=par(i)+edif;
    par2(i)=par(i)-edif;
    grad(i)=(functionhandle(par1)-functionhandle(par2))/(2*edif);
end