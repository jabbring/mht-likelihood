# The Likelihood of Mixed Hitting Times

This repository contains MATLAB code for replicating the numerical results in [Abbring and Salimans (2021)](https://arxiv.org/abs/1905.03463). This covers the computation of the likelihood of the [mixed hitting-time model](http://jaap.abbring.org/images/pdf/ecta7312.pdf), maximum likelihood estimation of parametric versions of this model, and an application to the analysis of [Kennan's (1985) strike data](https://www.ssc.wisc.edu/~jkennan/research/JEM85.pdf).

## Contents
To replicate [Abbring and Salimans (2021)](https://arxiv.org/abs/1905.03463), run the MATLAB script `replicate.m`. 

### Scripts for replication subtasks
These scripts are called by `replicate.m`. Adapt them to apply the code in other settings.
- `figure1.m`
- `figure2.m`
- `figure3.m`
- `simulateEstimate.m`
- `table1.m`
- `table1robustness.m`

### Functions for maximum likelihood estimation
- `mhtMaxLikelihood.m` - general case (based on Laplace inversion)
- `mhtMaxLikBM.m` - Gaussian special case

### Functions for Minus Log Likelihood calculation functions
- `minusLoglikelihood.m` - function that returns minus the log likelihood (calculated by Laplace transform inversion)
- `minusLoglikBM.m` - function that returns minus the log likelihood (calculated using explicit expressions for the Gaussian case)

### Specification functions
Each parametric specification of the model is coded up in a file `<heter><shocks>.m` as a function `<heter><shocks>`  that calculates the Laplace transform of .... 
for a model with unobserved heterogeneity specification `<heter>` and shock specification `<shocks>`:
- `pointpoint.m` - Discrete heterogeneity and discrete shocks at Poisson times
- `pointgamma.m` - Discrete heterogeneity and gamma shocks at Poisson times
- `gammapoint.m` - Gamma heterogeneity and discrete shocks at Poisson times
- `gammagamma.m` - Gamma heterogeneity and gamma shocks at Poisson times

Users can extend the set of specifications by adding different functions `<heter><shocks>`.

### Auxiliary functions
- `igausscdf.m`n - inverse Gaussian cdf
- `igausspdf.m` - inverse Gaussian pdf
- `randraw.m` -  random draws from ...

### Data
- `strdur.asc` - Fixed format text file with [Kennan's (1985) strike data](https://www.ssc.wisc.edu/~jkennan/research/JEM85.pdf) (source: [Cameron and Trivedi’s, 2005, data sets page](http://cameron.econ.ucdavis.edu/mmabook/mmadata.html)).


## References
- Abbring, Jaap H., and Tim Salimans (2021), “[The likelihood of mixed hitting times](https://arxiv.org/abs/1905.03463)”, *Journal of Econometrics*, forthcoming. arXiv:1905.03463 \[econ.EM\].
- Cameron, A. Colin, and Pravin K. Trivedi (2005), *[Microeconometrics: Methods and Applications](http://cameron.econ.ucdavis.edu/mmabook/mma.html)*, Cambridge: Cambridge University Press.
- Kennan, John (1985), "[The duration of contract strikes in U.S. manufacturing](https://www.ssc.wisc.edu/~jkennan/research/JEM85.pdf)", *Journal of Econometrics*, 28, 5–28.
