function pargrad = mhtgrad(mhtfun,par,y,cens,varargin)

[probs,grad]=numinvlap(mhtfun,par,y,cens,varargin{:});
pargrad=-grad'*(1./probs);