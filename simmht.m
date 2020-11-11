function y=simmht(n,x,mu,var,unobstype,v,p,shocktype,lambda,nu,beta,ymax)
%/////////////////////////////////////////////////////////////
% simulate durations from a MHT model
%/////////////////////////////////////////////////////////////

% set-up
y=zeros(n,1);
if length(ymax)==1
    ymax=ymax*ones(n,1);
end
gamma_unobs=isequal(unobstype,'gamma');
gamma_shock=isequal(shocktype,'gamma');
if gamma_unobs
    omega=v;
    tau=p;
    unobs=randraw('gamma',[0, 1/omega, tau],n);
else
    % unobserved heterogeneity
    unobs=v(1+sum(repmat(rand(n,1),1,length(p))>repmat(cumsum(p)',n,1),2));
end
if gamma_shock
    rho=nu(2);
    nu=nu(1);
else
    nu=abs(nu);
    lp=lambda/sum(lambda);
end
l=sum(lambda);
if ~isempty(x)
    unobs=unobs.*exp(x*beta);
end
par=[unobs/mu unobs.^2/var]';

% brownian motion + shocks
for i=1:n

    % start with brownian motion
    y(i)=randraw('ig',par(:,i),1);

    % integrate poisson process
    etime=y(i);
    while etime>0 && y(i)<ymax(i)
        nrshocks=randraw('po',l*etime,1);
        if nrshocks>0
            if gamma_shock
                edist=sum(randraw('gamma',[0, 1/nu, rho],nrshocks));
            else
                edist=sum(nu(1+sum(repmat(rand(nrshocks,1),1,length(lp))>repmat(cumsum(lp)',nrshocks,1),2)));
            end
            etime=randraw('ig',[edist/mu; edist^2/var],1);
            y(i)=y(i)+etime;
        else
            etime=0;
        end
    end
    if y(i)>ymax(i)
        y(i)=ymax(i);
    end

end
