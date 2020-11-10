# Matlab Code for "The Likelihood of Mixed Hitting Times"
## Jaap Abbring and Tim Salimans

This repository contains MATLAB code for computing the likelihood of a mixed hitting-time model that specifies durations as the first time a latent Lévy process crosses a heterogeneous threshold. This likelihood is not generally known in closed form, but its Laplace transform is. Our approach to its computation relies on numerical methods for inverting Laplace transforms that exploit special properties of the first passage times of Lévy processes. We use our method to implement a maximum likelihood estimator of the mixed hitting-time model. We illustrate the application of this estimator with an analysis of Kennan’s (1985) strike data.


