{smcl}
{cmd:help polbunch}{right: ()}
{hline}

{title:Title}

{p2colset 5 14 16 2}{...}
{p2col:{cmd:polbunch} {hline 2}}Theoretically consistent, model-based polynomial bunching estimation{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 13 2}
{cmd:polbunch} [{it:freqvar}] {it:zvar} {ifin}{cmd:,} {opt cut:off(#)} [{it:options}]

{pstd}
When one variable is specified, {it:zvar} is interpreted as individual-level earnings or log earnings, and {opt bw(#)} is required. When two
variables are specified, {it:freqvar} is interpreted as bin counts and {it:zvar} as bin midpoints; in that case the bandwidth is inferred
from the bin spacing and {opt bw()} may not be specified.


{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt cut:off(#)}}required; specifies the kink point or bunching point{p_end}
{synopt:{opt bw(#)}}bin width; required with individual-level data and not allowed with pre-binned data{p_end}
{synopt:{opt pol:ynomial(#)}}degree of the counterfactual polynomial; default is {cmd:polynomial(7)}{p_end}
{synopt:{opt lim:its(numlist)}}two integers specifying the number of excluded bins below and above the cutoff; default is {cmd:limits(1 0)}{p_end}
{synopt:{opt est:imator(#)}}estimator to use; default is {cmd:estimator(3)}{p_end}
{synopt:{opt log}}specifies that the running variable is in logs{p_end}
{synopt:{opt t0(#)}}linear tax rate below the cutoff{p_end}
{synopt:{opt t1(#)}}linear tax rate above the cutoff{p_end}

{syntab:Reporting and transformations}
{synopt:{opt notransform}}report raw estimating-equation coefficients rather than transformed bunching parameters{p_end}
{synopt:{opt constant}}use the constant-density approximation when transforming bunching into a response and elasticity{p_end}
{synopt:{opt Bmodel}}for estimator 2, report excess mass using the model-implied bunching mass rather than the observed excluded-bin mass{p_end}

{syntab:Estimation controls}
{synopt:{opt positive}}restrict the structural shift parameter (and therefore the elasticity) to be positive; by default shifts are restricted to be greater than -1{p_end}
{synopt:{opt nonormalize}}estimate on the original running-variable scale instead of normalizing bin centers by the cutoff and bin width{p_end}
{synopt:{opt nozero}}do not fill empty bins with zero counts; by default empty bins inside the observed support are included{p_end}
{synopt:{opt nodrop}}do not drop endpoint bins that appear to be cut by sample selection{p_end}
{synopt:{opt norankred}}do not reduce the polynomial degree when the unrestricted regression is rank-deficient; by default {cmd:polbunch} reduces the degree automatically{p_end}

{syntab:Inference}
{synopt:{opt vce(string)}}variance estimator; {cmd:analytic} (default), {cmd:bootstrap}, {cmd:bayes}, or {cmd:none}{p_end}
{synopt:{opt boot:reps(#)}}number of bootstrap repetitions when {cmd:vce(bootstrap)} or {cmd:vce(bayes)}; default is {cmd:bootreps(500)}; must be at least 2{p_end}
{synopt:{opt nodots}}suppress bootstrap progress dots{p_end}

{syntab:Model-restriction test}
{synopt:{opt test(string)}}which restriction test to report; {cmd:all} (default for estimators 2–3), {cmd:wald} (default for estimators 1 and 4), {cmd:hausman}, {cmd:minimumdistance}, or {cmd:none}{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:polbunch} estimates bunching at a concave kink using polynomial approximations to the counterfactual density. The command supports individual-level data, which it first collapses into bins,
and already binned data, where the first variable is the bin count and the second variable is the bin midpoint.

{pstd}
With individual-level data, {cmd:polbunch} forms bins of width {opt bw(#)} around the cutoff. Values exactly at the cutoff are assigned to the lower bin. Empty bins inside the observed support
are filled with zero counts unless {opt nozero} is specified. This makes the individual-level and pre-binned workflows comparable when they use the same support. If the support of the
individual-level data is much wider than the support of a pre-binned input, the estimates may differ because the polynomial is fit over a different range of bins; use {it:if} or {it:in}
to impose the desired estimation window.

{pstd}
By default, {cmd:polbunch} normalizes the estimation scale to

{p 12 12 2}
{it:z_est} = ({it:z_orig} - {it:cutoff})/{it:bw}.

{pstd}
This normalization improves numerical conditioning. Specify {opt nonormalize} to estimate on the original scale.

{pstd}
The option {opt limits(L H)} defines the excluded bunching region. If the cutoff lies inside a bin, that cutoff-crossing bin is included in the excluded region together with {it:L} bins below
and {it:H} bins above. If the cutoff lies exactly on a bin edge, there is no cutoff-crossing bin; the command excludes {it:L} bins below and {it:H} bins above the cutoff. Non-excluded control
bins are classified as left or right according to the edges of the excluded region.


{marker estimators}{...}
{title:Estimators}

{pstd}
{cmd:polbunch} implements five estimators. Estimators 0–3 use a unified polynomial/profile framework. Estimator 4 implements a separate Saez-style three-region trapezoid estimator.

{phang}
{cmd:estimator(0)} estimates an unrestricted model with separate left- and right-side polynomials and a free bunching mass. This estimator is useful for diagnostics and for testing the restrictions imposed by estimators 1–3.

{phang}
{cmd:estimator(1)} is the naive bunching estimator that does not correct for the distortions above the threshold. It sets the right-side counterfactual polynomial equal to the left-side polynomial. This estimator is generally biased under the isoelastic labor supply model.

{phang}
{cmd:estimator(2)} implements a Chetty et al. style adjustment. The right-side polynomial is proportional to the left-side polynomial and the bunching mass is tied to the implied missing area under the left-side counterfactual.

{phang}
{cmd:estimator(3)} is the default model-restricted polynomial bunching estimator. It imposes the density transformation implied by an isoelastic labor supply model. For level earnings, the right-side restriction is proportional in earnings;
for log earnings, it is additive in log earnings. The structural shift parameter is allowed to be greater than -1 by default; specify {opt positive} to require a positive shift.

{phang}
{cmd:estimator(4)} implements a Saez-style trapezoid approximation. It estimates a left reference height, a right reference height, and the bunching mass directly. This estimator does not use the nonlinear profile restrictions and does
not run the model-restriction test. With {opt notransform}, it reports {cmd:h0:_cons}, {cmd:h1:_cons}, and {cmd:bunching:B}. With the default transformation, it reports the reference heights, the number of bunchers, excess mass, and, when the trapezoid response equation can be solved, the shift, marginal response, and elasticity.

{pstd}
For estimators 2 and 3, {cmd:polbunch} uses a profile implementation of the corresponding nonlinear least-squares estimator. For any candidate value of the structural shift parameter, the remaining polynomial coefficients enter the model linearly
and are therefore estimated by least squares. The command then minimizes the resulting profiled sum of squared residuals over the structural parameter only. This recovers the same minimum and the same parameter estimates as the full nonlinear
least-squares problem, but avoids repeatedly optimizing over all polynomial coefficients jointly. In practice, the profile formulation is faster, less sensitive to starting values, and more numerically stable, especially for higher-order polynomials
or when the restricted model is tightly parameterized.
{p_end}

{marker transform}{...}
{title:Transformed parameters}

{pstd}
Unless {opt notransform} is specified, {cmd:polbunch} transforms the raw estimating-equation parameters into economically interpretable quantities. The transformed output may include the estimated counterfactual density under the low-tax regime {cmd:h0}, the density under the high-tax regime {cmd:h1}, the number of bunchers, excess mass, the proportional shift, the response of the marginal buncher, and, if {opt t0()} and {opt t1()} are specified, the elasticity.

{pstd}
For estimator 3, the response and elasticity are based on the model-implied density transformation. If {opt log} is specified, the response is a displacement in log earnings and the elasticity is the log response divided by the log net-of-tax-rate change. Without {opt log}, the level response is converted to a proportional shift before computing the elasticity.

{pstd}
For estimator 4, the Saez transformation first computes excess mass using the average of the two trapezoid endpoints. It then attempts to invert the trapezoid response equation. If the equation has no real positive solution, {cmd:polbunch} reports the reference heights, number of bunchers, and excess mass, but omits shift, marginal response, and elasticity and displays a note.

{pstd}
The {opt constant} option uses a constant-density approximation when converting bunching to a response and elasticity. This approximation can be useful for comparison with older procedures but may be biased when the density changes substantially over the response region.


{marker tests}{...}
{title:Model-restriction tests}

{pstd}
For estimators 1–4, {cmd:polbunch} tests the restrictions implied by the selected estimator, unless {cmd:test(none)} is specified or {cmd:vce(none)} is used. The test type is controlled by {opt test(string)}.

{phang}
{cmd:test(all)} (default for estimators 2 and 3) runs all applicable tests for the selected estimator and reports those that succeed. For estimators 2 and 3 this means minimum-distance and Hausman; for estimators 1 and 4 it reduces to the Wald test only.

{phang}
{cmd:test(hausman)} is a Hausman-type test that compares the restricted and unrestricted estimates. It uses the analytic variance of the difference between the two estimates when {cmd:vce(analytic)} is in effect, and the bootstrap covariance of the difference otherwise. Available for estimators 2 and 3.

{phang}
{cmd:test(wald)} (default for estimators 1 and 4) is a Wald test of the linear or nonlinear restrictions imposed by the selected estimator against the unrestricted estimator-0 estimates. For estimators 1 and 4 this reduces to a standard linear restriction test. When specified explicitly for estimators 2 or 3, {cmd:polbunch} prints a note recommending minimum-distance or Hausman instead, as the Wald statistic is a conditional shape diagnostic rather than a formal overall specification test.

{phang}
{cmd:test(minimumdistance)} is a minimum-distance test. It minimizes the Wald criterion over the structural parameter delta and compares the resulting minimum distance statistic to a chi-squared distribution. Available for estimators 2 and 3.

{phang}
{cmd:test(none)} suppresses all model-restriction testing.

{pstd}
All tests are reported as chi-squared statistics with associated p-values. They should be interpreted jointly as a test of the estimator's structural restrictions and the polynomial approximation used for the counterfactual density.


{marker inference}{...}
{title:Inference}

{pstd}
The variance estimator is controlled by {opt vce(string)}.

{phang}
{cmd:vce(analytic)} (the default) computes analytic standard errors using the collapsed-data delta-method variance correction. This correction reproduces the variance that would be obtained from a regression in the stacked data with one row per individual per bin and clustered standard errors,
without having to construct the expanded data. It is fast and can be used with pre-binned data.

{phang}
{cmd:vce(bootstrap)} performs the binned bootstrap. Instead of resampling individual observations directly, the command draws bin counts from Gamma distributions scaled to preserve the total sample size. This mimics the classical bootstrap while remaining feasible when only binned data are available. The number of repetitions is set by {opt bootreps(#)}; the default is 500 and the minimum is 2.

{phang}
{cmd:vce(bayes)} performs the Bayesian bootstrap using Dirichlet sampling of bin probabilities. Each repetition draws a new set of bin weights from the posterior implied by a Dirichlet prior and updates bin counts accordingly. The number of repetitions is set by {opt bootreps(#)}.

{phang}
{cmd:vce(none)} suppresses all internal variance estimation. Point estimates are computed but no standard errors or covariance matrix are stored. This is useful in Monte Carlo exercises where only point estimates are needed, or when using Stata's {cmd:bootstrap} prefix around {cmd:polbunch}.

{pstd}
The analytic and bootstrap procedures treat individuals within a bin as identical from the point of view of the estimator and assume the bandwidth and binning scheme are fixed. If inference must account for clustering at a higher level, use a resampling procedure outside {cmd:polbunch}.


{marker options_detail}{...}
{title:Details on selected options}

{phang}
{opt log} tells {cmd:polbunch} that the running variable is already in logs. It changes the model-implied response mapping and the elasticity transformation. It does not log-transform the variable for the user.

{phang}
{opt t0(#)} and {opt t1(#)} specify the tax rates below and above the cutoff. If one is specified, the other must also be specified. The rates must differ, and the current implementation requires {cmd:t1()>t0()} (a convex kink).

{phang}
{opt positive} restricts the estimator-3 structural shift to be positive. Without {opt positive}, estimator 3 allows negative shifts as long as {cmd:1 + delta > 0}.

{phang}
{opt nonormalize} leaves the running variable on its original scale during estimation. The default normalization usually improves numerical conditioning and does not change the population estimand.

{phang}
{opt nozero} affects only individual-level input. By default, after collapsing individual observations into bins, empty bins inside the observed support are retained as zero-count bins. {opt nozero} excludes such bins from estimation.

{phang}
{opt Bmodel} affects transformed results for estimator 2. It requests the model-implied bunching mass rather than the observed excluded-bin mass when constructing transformed bunching quantities.

{phang}
{opt norankred} suppresses the automatic polynomial degree reduction that {cmd:polbunch} performs when the unrestricted two-sided polynomial regression is rank-deficient. By default, the degree is reduced one step at a time until identification is restored; a note is displayed.

{phang}
{opt test(string)} selects the restriction test. The default is {cmd:all} for estimators 2 and 3 and {cmd:wald} for estimators 1 and 4. Specifying {cmd:test(none)} skips all testing. Note that testing is also suppressed when {cmd:vce(none)} is in effect.


{marker examples}{...}
{title:Examples}

{pstd}
Generate simulated data from the companion data-generating command:{p_end}

{phang2}{cmd:. polbunchgendata, obs(10000) t0(0.2) t1(0.6) el(0.4) cutoff(1)}{p_end}

{pstd}
Estimate bunching using the default estimator and a correctly specified first-degree polynomial with bandwidth 0.01:{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(1) t0(0.2) t1(0.6)}{p_end}

{pstd}
Use the same estimator without transforming to economic bunching parameters:{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(1) estimator(3) notransform}{p_end}

{pstd}
Estimate on the original running-variable scale rather than the normalized scale:{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(1) nonormalize}{p_end}

{pstd}
Compare with the Chetty-style adjustment estimator:{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(1) estimator(2)}{p_end}

{pstd}
Compare with the naive no-adjustment estimator:{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(1) estimator(1)}{p_end}

{pstd}
Compare with the Saez trapezoid estimator, restricting to a small region around the cutoff:{p_end}

{phang2}{cmd:. polbunch z if inrange(z,0.9,1.1), cutoff(1) bw(0.01) estimator(4) t0(0.2) t1(0.6)}{p_end}

{pstd}
Use the binned bootstrap for inference with 200 repetitions:{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(1) vce(bootstrap) bootreps(200)}{p_end}

{pstd}
Suppress internal variance estimation, for example when using Stata's bootstrap prefix:{p_end}

{phang2}{cmd:. bootstrap, reps(200): polbunch z, cutoff(1) bw(0.01) polynomial(1) vce(none)}{p_end}

{pstd}
Request a minimum-distance restriction test instead of the default Hausman test:{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(1) test(minimumdistance)}{p_end}

{pstd}
Suppress all model-restriction testing:{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(1) test(none)}{p_end}

{pstd}
Collapse to binned data and use polbunch with bin counts {cmd:freq} and bin midpoints {cmd:zmid}:{p_end}

{phang2}{cmd:. gen bin = ceil((z-1)/.01)*.01 + 1 - .005}{p_end}
{phang2}{cmd:. collapse (count) freq=z, by(bin)}{p_end}
{phang2}{cmd:. polbunch freq bin, cutoff(1) pol(1)}{p_end}


{marker saved_results}{...}
{title:Stored results}

{pstd}
{cmd:polbunch} stores the following in {cmd:e()}:

{synoptset 26 tabbed}{...}
{p2col 5 26 30 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations represented by the estimation sample{p_end}
{synopt:{cmd:e(estimator)}}estimator number{p_end}
{synopt:{cmd:e(polynomial)}}degree of polynomial used{p_end}
{synopt:{cmd:e(cutoff_orig)}}cutoff on the original running-variable scale{p_end}
{synopt:{cmd:e(cutoff_est)}}cutoff on the estimation scale (after normalization){p_end}
{synopt:{cmd:e(bw)}}bin width on the original scale{p_end}
{synopt:{cmd:e(bw_orig)}}bin width on the original scale{p_end}
{synopt:{cmd:e(bw_est)}}bin width on the estimation scale{p_end}
{synopt:{cmd:e(xscale)}}normalization scale factor ({it:bw_orig} / {it:bw_est}){p_end}
{synopt:{cmd:e(lower_limit)}}lower edge of the excluded region on the original scale{p_end}
{synopt:{cmd:e(upper_limit)}}upper edge of the excluded region on the original scale{p_end}
{synopt:{cmd:e(zL_excl_est)}}lower edge of the excluded region on the estimation scale{p_end}
{synopt:{cmd:e(zH_excl_est)}}upper edge of the excluded region on the estimation scale{p_end}
{synopt:{cmd:e(dL)}}mean bin midpoint in the left reference region{p_end}
{synopt:{cmd:e(dR)}}mean bin midpoint in the right reference region{p_end}
{synopt:{cmd:e(zlo)}}bottom of the lowest bin in the estimation window (original scale){p_end}
{synopt:{cmd:e(zhi)}}top of the highest bin in the estimation window (original scale){p_end}
{synopt:{cmd:e(log)}}1 if {opt log} was specified, 0 otherwise{p_end}
{synopt:{cmd:e(chi2_wald)}}Wald test chi-squared statistic, when computed{p_end}
{synopt:{cmd:e(p_wald)}}Wald test p-value, when computed{p_end}
{synopt:{cmd:e(df_wald)}}Wald test degrees of freedom, when computed{p_end}
{synopt:{cmd:e(chi2_minimumdistance)}}minimum-distance test chi-squared statistic, when computed{p_end}
{synopt:{cmd:e(p_minimumdistance)}}minimum-distance test p-value, when computed{p_end}
{synopt:{cmd:e(df_minimumdistance)}}minimum-distance test degrees of freedom, when computed{p_end}
{synopt:{cmd:e(chi2_hausman)}}Hausman test chi-squared statistic, when computed{p_end}
{synopt:{cmd:e(p_hausman)}}Hausman test p-value, when computed{p_end}
{synopt:{cmd:e(df_hausman)}}Hausman test degrees of freedom, when computed{p_end}
{synopt:{cmd:e(delta_md)}}structural shift estimate from the minimum-distance test, when computed{p_end}

{synoptset 26 tabbed}{...}
{p2col 5 26 30 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:polbunch}{p_end}
{synopt:{cmd:e(cmdname)}}{cmd:polbunch}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}dependent variable name in the estimation output{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(binname)}}name of the running-variable/bin variable{p_end}
{synopt:{cmd:e(normalize)}}normalization flag ({cmd:nonormalize} if specified, otherwise empty){p_end}
{synopt:{cmd:e(transform)}}transformation flag ({cmd:notransform} if specified, otherwise empty){p_end}
{synopt:{cmd:e(vcetype)}}variance-estimator label, when set{p_end}
{synopt:{cmd:e(properties)}}usually {cmd:b V} when a variance matrix is posted{p_end}

{synoptset 26 tabbed}{...}
{p2col 5 26 30 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector. With default transformation this contains density and bunching parameters; with {opt notransform} it contains raw estimating-equation coefficients{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of {cmd:e(b)}, when computed{p_end}
{synopt:{cmd:e(table)}}table of binned frequencies and variables used for plotting and diagnostics{p_end}
{synopt:{cmd:e(G)}}delta-method Jacobian for transformed parameters, when available{p_end}

{synoptset 26 tabbed}{...}
{p2col 5 26 30 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks the estimation sample{p_end}


{marker references}{...}
{title:References}

{phang}
Andresen, Martin E. (2026). "A better polynomial bunching estimator", working paper.

{phang}
Kleven, Henrik Jacobsen (2016). "Bunching", {it:Annual Review of Economics}.

{phang}
Saez, Emmanuel (2010). "Do Taxpayers Bunch at Kink Points?", {it:American Economic Journal: Economic Policy}.

{phang}
Chetty, Raj, John N. Friedman, Tore Olsen, and Luigi Pistaferri (2011). "Adjustment Costs, Firm Responses, and Micro vs. Macro Labor Supply Elasticities: Evidence from Danish Tax Records", {it:Quarterly Journal of Economics}.


{title:Suggested citation}

{pstd}
Andresen, Martin E. (2026). "POLBUNCH: Stata module for the polynomial bunching estimator." This version VERSION_DATE.{p_end}

{pstd}
Check your installed version date with:{p_end}

{phang2}{cmd:. which polbunch}{p_end}


{marker author}{...}
{title:Author}

{pstd}Martin Eckhoff Andresen{p_end}
{pstd}University of Oslo{p_end}
{pstd}Department of Economics{p_end}
{pstd}Oslo, Norway{p_end}
{pstd}martin.eckhoff.andresen@gmail.com{p_end}


{marker also_see}{...}
{title:Also see}

{p 4 14 2}
Development version: net install polbunch, from("https://raw.githubusercontent.com/martin-andresen/polbunch/master/"){p_end}

{p 7 14 2}
Help: {helpb polbunchplot}, {helpb polbunchgendata}, {helpb polbunchsim}{p_end}
