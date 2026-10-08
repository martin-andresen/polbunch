# polbunch

**Bunching estimation consistent with the model it assumes.**

Conventional polynomial bunching estimators (Saez 2010; Chetty et al. 2011) ignore or only approximately handle how the iso-elastic response to a kink reshapes the density *above* the cutoff, leading to bias. `polbunch` implements an estimator that is consistent and efficient under the iso-elastic model as well as the existing estimators from the literature, requiring only binned data.

<img src="polbunchex.png" alt="Bunching plot" width="600"/>

## Features

- **Consistent and efficient estimator.** The default `estimator(3)` imposes the density transformation implied by the iso-elastic model and estimates elasticity, excess mass and counterfactual jointly, exploiting all the restrictions of the assumed model.
- **Existing estimators.** Saez-style (`estimator(4)`), Chetty-style (`estimator(5)`, and the naive estimator (`estimator(1)`) run through the same interface, on the same bins, so results are directly comparable. 
- **Contrast results** `polbunch_contrast` tests whether two specifications really give different elasticities, accounting for their dependence.
- **Bias characterization.** `polbunchbias` computes the asymptotic bias of a conventional estimator for a given true elasticity, tax change and counterfactual density, without simulation. With already biased inputs, iterate the procedure for a fixed point.
- **Analytic standard errors.** Closed-form `vce(analytic)` works with pre-binned data, so no micro data are needed. Binned bootstrap variants (multinomial, residual, wild, Bayesian) also available. Analytic clustered standard errors when micro-data is available or the co-visitation matrix exist alongside the histogram for binned data.
- **Specification tests.** Compare the restrictions your estimator imposes against an unrestricted two-sided fit (minimum-distance, omnibus, Hausman). `estat gof` gives reference-region goodness of fit, including separate fit below and above the cutoff.
- **Permutation inference.** `polbunch_permute` runs a placebo-cutoff permutation test on the studentized excess mass.
- **Plotting tools.** `polbunchplot` draws the observed and counterfactual densities for one or multiple polbunch estimateds, or a rootogram to show where the fit fails.

## Quick start

```stata
* Generate iso-elastic data
polbunchgendata z

* Estimate with the efficient estimator and analytic standard errors
polbunch z, cutoff(1) bw(0.01) polynomial(5) estimator(3) t0(0.2) t1(0.6) 

* Plot observed vs. counterfactual density
polbunchplot

* How wrong would the Chetty-style estimator be here?
polbunch z, cutoff(1) bw(0.01) polynomial(5) estimator(5) t0(0.3) t1(0.4) contrast vce(bootstrap)

* Goodness of fit and rootogram
estat gof
polbunchplot, rootogram
```


## Install

```stata
net install polbunch, from("https://raw.githubusercontent.com/martin-andresen/polbunch/master/") replace
```

Then see `help polbunch`, `help polbunchbias`, `help polbunch_permute`, `help polbunch_contrast`, `help polbunchplot`.

## Citation

Andresen, M. E. (2026) "An Elasticity From a Single Histogram? Bias, Inference and Specification Tests for Polynomial Bunching Estimators." Working paper.

## Status

In development; bugs are possible. Report them through GitHub issues or to martin.eckhoff.andresen@gmail.com.
