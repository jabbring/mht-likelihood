function cdf=weibullcdf(t,a,d);
    cdf = 1-exp(-(t.^d).*a);
end