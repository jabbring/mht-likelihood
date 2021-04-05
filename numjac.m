function jac=numjac(functionhandle,par)
edif=1e-6;
jac=[];
for j=1:numel(par)
    par1=par;
    par2=par;
    par1(j)=par(j)+edif;
    par2(j)=par(j)-edif;
    jac=[jac 0.5*(functionhandle(par1)-functionhandle(par2))/edif];
end