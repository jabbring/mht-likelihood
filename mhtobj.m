function [nllh, pargrad] = mhtobj(mhtfun,par,y,cens,varargin)

[probs,grad]=numinvlap(mhtfun,par,y,cens,varargin{:});
nllh=(~cens)'*log(y)-sum(log(probs));
if nargout==2
    pargrad=-grad'*(1./probs);
end