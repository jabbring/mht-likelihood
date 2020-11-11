function [pd, grad]=igausspdf(y,mu,lambda)
pd=sqrt(lambda./(2*pi*(y.^3))).*exp(-lambda.*((y-mu).^2)./((2*mu.^2).*y));
if nargout==2
    grad=[pd.*(2*lambda.*(y-mu)./(2*mu.^2.*y)+4*mu.*y.*lambda.*(y-mu).^2./(2*mu.^2.*y).^2) ...
        pd.*(0.5./lambda-(y-mu).^2./(2*mu.^2.*y))];
end