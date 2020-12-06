function [nllh, pargrad] = mhtobj(mhtfun,par,y,M,cens,varargin)

[probs,grad]=numinvlap2(mhtfun,par,y,M,cens,varargin{:});
nllh=(~cens)'*log(y)-sum(log(probs));
if nargout==2
    pargrad=-grad'*(1./probs);
end