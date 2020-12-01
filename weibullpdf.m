function [pdf,grad]=weibullpdf(t,a,d);
    pdf = d.*(t.^(d-1)).*a.*exp(-(t.^d).*a);
    if nargout==2
        grad=[((1./a)-t.^d).*pdf ((1./d)+log(t).*(1-(t.^d).*a)).*pdf];
    end
end