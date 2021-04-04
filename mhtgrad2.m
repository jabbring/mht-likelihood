function pargrad = mhtgrad(mhtfun,par,y,M,cens,varargin)

[probs,grad]=numinvlap2(mhtfun,par,y,M,cens,varargin{:});
pargrad=-grad'*(1./probs);