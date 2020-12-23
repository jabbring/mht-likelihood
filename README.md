# The Likelihood of Mixed Hitting Times: Replication Package

This repository contains MATLAB code for replicating the numerical results in [Abbring and Salimans (2021)](https://arxiv.org/abs/1905.03463). This covers the computation of the likelihood of the [mixed hitting-time model](http://jaap.abbring.org/images/pdf/ecta7312.pdf), maximum likelihood estimation of parametric versions of this model, and an application to the analysis of [Kennan's (1985) strike data](https://www.ssc.wisc.edu/~jkennan/research/JEM85.pdf). 

The results in [Abbring and Salimans (2021)](https://arxiv.org/abs/1905.03463) were obtained using [v1.1.0](https://github.com/jabbring/mht-likelihood/releases/tag/v1.1.0) of this package, by running `make` on a MacBook Pro (2018, 15inch, 2.9GHz 6-Core Intel Core i9, 32 GB 2400 MHz DDR4) with macOS 10.15.7, MATLAB R2020b (with Optimization Toolbox), and TeX Live 2018. This executes the MATLAB scripts

- [figure1.m](https://github.com/jabbring/mht-likelihood/blob/master/figure1.m) -  replicates Figure 1
- [figure2.m](https://github.com/jabbring/mht-likelihood/blob/master/figure2.m) -  replicates Figure 2
- [figure3.m](https://github.com/jabbring/mht-likelihood/blob/master/figure3.m) -  replicates Figure 3
- [table1.m](https://github.com/jabbring/mht-likelihood/blob/master/table1.m) - replicates Table 1
- [table1lowM.m](https://github.com/jabbring/mht-likelihood/blob/master/table1lowM.m) - recalculates Table 1 with *M=15*
- [table1mig.m](https://github.com/jabbring/mht-likelihood/blob/master/table1mig.m) - recalculates Columns I-V of table one using the exact likelihood for the Gaussian case
- [figure4.m](https://github.com/jabbring/mht-likelihood/blob/master/figure4.m) - replicates Figure 4
- [checkgradient.m](https://github.com/jabbring/mht-likelihood/blob/master/checkgradient.m) - checks the analytic gradients against numerical ones

and displays the resulting figures and tables in a `pdf` file [replication.pdf](https://github.com/jabbring/mht-likelihood/blob/master/replication.pdf) (note that TeX Live is only used in the last step, which requires pdfTeX to process the `tex` files written by the MATLAB scripts). As an extra check, the scripts for computing Figures 1-4 can be adapted, by changing `dispplot=false` into `dispplot=true`, to plot these figures directly. The repository also contains a script

- [simtest.m](https://github.com/jabbring/mht-likelihood/blob/master/simtest.m) - simulates data and estimates the model on these data

Users can adapt these scripts to apply the procedures they call in other contexts. In particular, `table1.m` shows how to use the package to estimate the MHT model with jumps. 

The scripts require a range of functions (1-5) and a data set (6).

### 1. Model specification

The procedures require parametric specifications `<heter>` of the unobserved heterogeneity specification and `<shocks>` of the jumps in the latent Lévy process. The calculation of the Laplace transform of the mixed survival function for each such specification is coded up as a function `<heter><shocks>` in a file `<heter><shocks>.m`:

- [pointpoint.m](https://github.com/jabbring/mht-likelihood/blob/master/pointpoint.m) - Discrete heterogeneity and discrete shocks at Poisson times
- [pointgamma.m](https://github.com/jabbring/mht-likelihood/blob/master/pointgamma.m) - Discrete heterogeneity and gamma shocks at Poisson times
- [gammapoint.m](https://github.com/jabbring/mht-likelihood/blob/master/gammapoint.m) - Gamma heterogeneity and discrete shocks at Poisson times
- [gammagamma.m](https://github.com/jabbring/mht-likelihood/blob/master/gammagamma.m)- Gamma heterogeneity and gamma shocks at Poisson times

Users can extend the set of specifications by adding different functions `<heter><shocks>`.

### 2. Probability densities and cumulate distributions

- [numinvlap.m](https://github.com/jabbring/mht-likelihood/blob/master/numinvlap.m) - MHT pdf (calculated by Laplace transform inversion)
- [numinvlap2.m](https://github.com/jabbring/mht-likelihood/blob/master/numinvlap2.m) - idem, but with input of design parameter `M`
- [igausscdf.m](https://github.com/jabbring/mht-likelihood/blob/master/igausscdf.m) - inverse Gaussian cdf
- [igausspdf.m](https://github.com/jabbring/mht-likelihood/blob/master/igausspdf.m) - inverse Gaussian pdf
- [weibullcdf.m](https://github.com/jabbring/mht-likelihood/blob/master/weibullcdf.m) - Weibull MPH cdf
- [weibullpdf.m](https://github.com/jabbring/mht-likelihood/blob/master/weibullpdf.m) - Weibull MPH pdf

### 3. Likelihood calculation

- [mhtobj.m](https://github.com/jabbring/mht-likelihood/blob/master/mhtobj.m) - minus the log likelihood of the MHT model (calculated by Laplace transform inversion)
- [lhmigauss.m](https://github.com/jabbring/mht-likelihood/blob/master/lhmigauss.m) - likelihood of the MHT model (calculated using explicit expressions for the Gaussian case)
- [nllhmigauss.m](https://github.com/jabbring/mht-likelihood/blob/master/nllhmigauss.m) - minus the log likelihood of the MHT model (calculated using explicit expressions for the Gaussian case)
- [nllhmph.m](https://github.com/jabbring/mht-likelihood/blob/master/nllhmph.m) - minus the log likelihood of the Weibull MPH model

### 4. Maximum likelihood estimation

- [mhtmle.m](https://github.com/jabbring/mht-likelihood/blob/master/mhtmle.m) - ML estimation of the general MHT model (based on Laplace inversion)
- [migaussmle.m](https://github.com/jabbring/mht-likelihood/blob/master/migaussmle.m) - ML estimation of the Gaussian special case

### 5. Simulation

- [randraw.m](https://github.com/jabbring/mht-likelihood/blob/master/randraw.m) -  efficient random number generator (by Alex Bar Guy and Alexander Podgaetsky)
- [simmht.m](https://github.com/jabbring/mht-likelihood/blob/master/simmht.m) - simulates durations from mixed hitting time model

### 6. Data

- [strkdur.asc](https://github.com/jabbring/mht-likelihood/blob/master/strkdur.asc) - Fixed format text file with [Kennan's (1985) strike data](https://www.ssc.wisc.edu/~jkennan/research/JEM85.pdf) (source: [Cameron and Trivedi’s, 2005, data sets page](http://cameron.econ.ucdavis.edu/mmabook/mmadata.html)).

## References
- Abbring, Jaap H., and Tim Salimans (2021), “[The likelihood of mixed hitting times](https://arxiv.org/abs/1905.03463)”, *Journal of Econometrics*, forthcoming. arXiv:1905.03463 \[econ.EM\].
- Cameron, A. Colin, and Pravin K. Trivedi (2005), *[Microeconometrics: Methods and Applications](http://cameron.econ.ucdavis.edu/mmabook/mma.html)*, Cambridge: Cambridge University Press.
- Kennan, John (1985), "[The duration of contract strikes in U.S. manufacturing](https://www.ssc.wisc.edu/~jkennan/research/JEM85.pdf)", *Journal of Econometrics*, 28, 5–28.

We welcome the use of this software under an [MIT license](https://github.com/jabbring/mht-likelihood/blob/master/LICENSE).

&copy; 2020 Jaap H. Abbring and Tim Salimans
