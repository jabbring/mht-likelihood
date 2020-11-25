function [cd, grad]=igausscdf(y,mu,lambda)
cd=0.5*erfc(-sqrt(lambda./y).*(y./mu-1)/sqrt(2))+exp(2*lambda./mu).*0.5.*erfc(sqrt(lambda./y).*(y./mu+1)./sqrt(2));
if nargout==2
    grad=[-exp(-(lambda./y).*(y./mu-1).^2/2)./sqrt(pi).*(sqrt(lambda./y).*(y./mu.^2)/sqrt(2))-...
        exp(2*lambda./mu).*(lambda./mu.^2).*erfc(sqrt(lambda./y).*(y./mu+1)./sqrt(2))+...
        exp(2*lambda./mu).*exp(-(lambda./y).*(y./mu+1).^2./2).*(sqrt(lambda./y).*(y./mu.^2)./sqrt(2))./sqrt(pi) ...
        exp(-(lambda./y).*(y./mu-1).^2/2)./sqrt(pi).*(0.5*sqrt(y./lambda).*(1./y).*(y./mu-1)/sqrt(2))+...
        exp(2*lambda./mu).*(1./mu).*erfc(sqrt(lambda./y).*(y./mu+1)./sqrt(2))-...
        exp(2*lambda./mu).*exp(-(lambda./y).*(y./mu+1).^2./2).*(0.5*sqrt(y./lambda).*(1./y).*(y./mu+1)./sqrt(2))/sqrt(pi)];
end