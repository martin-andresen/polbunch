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
{synopt:{opt lim:its(numlist)}}two integers specifying the number of excluded bins below and above the cutoff; default is {cmd:limits(1 0)} L>0, H>=0.{p_end}
{synopt:{opt est:imator(#)}}estimator to use; default is {cmd:estimator(3)}{p_end}
{synopt:{opt log}}specifies that the running variable is in logs{p_end}
{synopt:{opt t0(#)}}linear tax rate below the cutoff{p_end}
{synopt:{opt t1(#)}}linear tax rate above the cutoff{p_end}

{syntab:Reporting and transformations}
{synopt:{opt notransform}}report raw estimating-equation coefficients rather than transformed bunching parameters{p_end}
{synopt:{opt constant}}use the constant-density approximation when transforming bunching into a response and elasticity{p_end}
{synopt:{opt Bmodel}}for estimator 2, report excess mass using the model-implied bunching mass rather than the observed excluded-bin mass{p_end}

{syntab:Estimation controls}
{synopt:{opt positive}}restrict the structural shift parameter (and therefore the elasticity) to be positive; by default shfifts are restricted to be >-1{p_end}
{synopt:{opt nonormalize}}estimate on the original running-variable scale instead of normalizing bin centers by the cutoff and bin width{p_end}
{synopt:{opt nozero}}do not fill empty bins with zero counts before estimation or during bootstrap - by default empty bins are included. {p_end}
{synopt:{opt nodrop}}do not drop endpoint bins that appear to be cut by sample selection{p_end}

{syntab:Inference and diagnostics}
{synopt:{opt boot:reps(#)}}number of binned bootstrap repetitions; the default {cmd:bootreps(1)} gives analytic standard errors, {cmd:bootreps(0)} suppresses internal variance estimation{p_end}
{synopt:{opt notest}}do not test model restrictions against the unrestricted polynomial model{p_end}
{synopt:{opt saveunres(name)}}store the unrestricted comparison model used for the restriction test{p_end}
{synopt:{opt nodots}}suppress bootstrap progress dots{p_end}
{synopt:{opt nosmallsample}}omit the small-sample adjustment in the analytic variance correction{p_end}
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
{cmd:polbunch} implements five estimators. Estimators 0--3 use a unified polynomial/profile framework. Estimator 4 implements a separate Saez-style three-region trapezoid estimator.

{phang}
{cmd:estimator(0)} estimates an unrestricted model with separate left- and right-side polynomials and a free bunching mass. This estimator is useful for diagnostics and for testing the restrictions imposed by estimators 1--3.

{phang}
{cmd:estimator(1)} is the naive bunching estimator that does not correct for the distortions above the threshold. It sets the right-side counterfactual polynomial equal to the left-side polynomial. This estimator is generally biased under the isoelastic labor supply model.

{phang}
{cmd:estimator(2)} implements a Chetty et al. style adjustment. The right-side polynomial is proportional to the left-side polynomial and the bunching mass is tied to the implied missing area under the left-side counterfactual.

{phang}
{cmd:estimator(3)} is the default model-restricted polynomial bunching estimator. It imposes the density transformation implied by an isoelastic labor supply model. For level earnings, the right-side restriction is proportional in earnings; 
for log earnings, it is additive in log earnings. The structural shift parameter is allowed to be greater than -1 by default; specify {opt positive} to require a positive shift.

{phang}
{cmd:estimator(4)} implements a Saez-style trapezoid approximation. It estimates a left reference height, a right reference height, and the bunching mass directly. This estimator does not use the nonlinear profile restrictions and does 
not run the model-restriction test. With {opt notransform}, it reports {cmd:h0:_cons}, {cmd:h1:_cons}, and {cmd:bunching:B}. With the default transformation, it reports the reference heights, the number of bunchers, excess mass, and, 
when the trapezoid response equation can be solved, the shift, marginal response, and elasticity.


{marker transform}{...}
{title:Transformed parameters}

{pstd}
Unless {opt notransform} is specified, {cmd:polbunch} transforms the raw estimating-equation parameters into economically interpretable quantities. The transformed output may include the estimated counterfactual density under the low-tax regime {cmd:h0}, the density under the high-tax regime {cmd:h1}, the number of bunchers, excess mass, the proportional shift, the response of the marginal buncher, and, if {opt t0()} and {opt t1()} are specified, the elasticity.

{pstd}
For estimator 3, the response and elasticity are based on the model-implied density transformation. If {opt log} is specified, the response is a displacement in log earnings and the elasticity is the log response divided by the log net-of-tax-rate change. Without {opt log}, the level response is converted to a proportional shift before computing the elasticity.

{pstd}
For estimator 4, the Saez transformation first computes excess mass using the average of the two trapezoid endpoints. It then attempts to invert the trapezoid response equation. If the equation has no real positive solution, {cmd:polbunch} reports the reference heights, number of bunchers, and excess mass, but omits shift, marginal response, and elasticity and displays a note. This can occur, for example, when the estimated bunching mass is too large relative to the endpoint heights and excluded-region width.

{pstd}
The {opt constant} option uses a constant-density approximation when converting bunching to a response and elasticity. This approximation can be useful for comparison with older procedures but may be biased when the density changes substantially over the response region.


{marker tests}{...}
{title:Model-restriction tests}

{pstd}
For estimators 1--3, {cmd:polbunch} additionally estimates the unrestricted estimator 0 and tests the restrictions implied by the selected estimator, unless {opt notest} is specified or internal variance estimation is suppressed. The test is reported as a 
chi-squared statistic with a p-value. The test should be interpreted jointly as a test of the estimator's structural restrictions and the polynomial approximation used for the counterfactual density.

{pstd}
Estimator 4 is a separate Saez-style comparison estimator and does not impose the polynomial profile restrictions tested for estimators 1--3.


{marker inference}{...}
{title:Inference}

{pstd}
The default {cmd:bootreps(1)} computes analytic standard errors using the collapsed-data variance correction. This correction reproduces the variance that would be obtained from an imaginary regression in the stacked data with one row per individual and bin and clustered standard errors,
without having to construct the expanded data. It can be used with binned data and is fast in simulations.

{pstd}
If {cmd:bootreps(}{it:B}{cmd:)} is specified with {it:B} > 1, {cmd:polbunch} performs the binned bootstrap. Instead of resampling individual observations directly, the command resamples bin counts from multinomial or binomial sampling probabilities implied by the observed binned distribution. This is designed to mimic the classical bootstrap while remaining feasible when only binned data are available.

{pstd}
Specify {cmd:bootreps(0)} to suppress internal variance estimation. This is useful in Monte Carlo exercises where only point estimates are needed, or when using Stata's {cmd:bootstrap} prefix around {cmd:polbunch}.

{pstd}
The analytic and binned-bootstrap procedures treat individuals within a bin as identical from the point of view of the estimator and assume the bandwidth and binning scheme are fixed. If inference must account for clustering at a higher level, use a resampling procedure outside {cmd:polbunch}.


{marker options_detail}{...}
{title:Details on selected options}

{phang}
{opt log} tells {cmd:polbunch} that the running variable is already in logs. It changes the model-implied response mapping and the elasticity transformation. It does not log-transform the variable for the user.

{phang}
{opt t0(#)} and {opt t1(#)} specify the tax rates below and above the cutoff. If one is specified, the other must also be specified. The rates must differ, and the current implementation is intended for convex kinks with {cmd:t1()>t0()}.

{phang}
{opt positive} restricts the estimator 3 structural shift to be positive. Without {opt positive}, estimator 3 allows negative shifts as long as {cmd:1 + delta > 0}.

{phang}
{opt nonormalize} leaves the running variable on its original scale during estimation. The default normalization usually improves numerical conditioning and should not change the population estimand.

{phang}
{opt nozero} affects only individual-level input. By default, after collapsing individual observations into bins, empty bins inside the observed support are retained as zero-count bins. {opt nozero} leaves such bins out of the estimation data.

{phang}
{opt Bmodel} affects transformed results for estimator 2. It requests the model-implied bunching mass rather than the observed excluded-bin mass when constructing transformed bunching quantities.

{phang}
{opt saveunres(name)} stores the unrestricted estimator 0 model used for the model-restriction test under the supplied name.


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
Use the binned bootstrap for inference:{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(1) bootreps(200)}{p_end}

{pstd}
Collapse to binned data and use polbunch with bin counts {cmd:freq} and bin midpoints {cmd:zmid}:{p_end}

{phang2}{cmd:. gen bin = ceil((z-1)/.01)*.01 + 1 - .005}{p_end}
{phang2}{cmd:. collapse (count) freq=z, by(bin)}{p_end}
{phang2}{cmd:. polbunch freq bin, cutoff(1) pol(1)}{p_end}

{pstd}
Suppress internal variance estimation, for example when using Stata's bootstrap prefix:{p_end}

{phang2}{cmd:. bootstrap, reps(200): polbunch z, cutoff(1) bw(0.01) polynomial(1) bootreps(0)}{p_end}


{marker saved_results}{...}
{title:Stored results}

{pstd}
{cmd:polbunch} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations represented by the estimation sample{p_end}
{synopt:{cmd:e(estimator)}}estimator number{p_end}
{synopt:{cmd:e(polynomial)}}degree of polynomial used{p_end}
{synopt:{cmd:e(cutoff)}}cutoff on the original scale{p_end}
{synopt:{cmd:e(bw)}}bin width on the original scale{p_end}
{synopt:{cmd:e(bandwidth)}}bin width on the original scale{p_end}
{synopt:{cmd:e(lower_limit)}}lower edge of the excluded region on the original scale{p_end}
{synopt:{cmd:e(upper_limit)}}upper edge of the excluded region on the original scale{p_end}
{synopt:{cmd:e(log)}}1 if {opt log} was specified, 0 otherwise{p_end}
{synopt:{cmd:e(chi2)}}chi-squared statistic for the model-restriction test, when computed{p_end}
{synopt:{cmd:e(p_mod)}}p-value for the model-restriction test, when computed{p_end}
{synopt:{cmd:e(df_mod)}}degrees of freedom for the model-restriction test, when computed{p_end}
{synopt:{cmd:e(hasresp)}}for estimator 4 after transformation, 1 if the Saez response equation was solved and 0 otherwise{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:polbunch}{p_end}
{synopt:{cmd:e(cmdname)}}{cmd:polbunch}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}dependent variable name in the estimation output{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(binname)}}name of the running-variable/bin variable{p_end}
{synopt:{cmd:e(normalize)}}normalization status{p_end}
{synopt:{cmd:e(vcetype)}}variance-estimator label, when set{p_end}
{synopt:{cmd:e(properties)}}usually {cmd:b V} when a variance matrix is posted{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector. With default transformation this contains density and bunching parameters; with {opt notransform} it contains raw estimating-equation coefficients{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of {cmd:e(b)}, when computed{p_end}
{synopt:{cmd:e(table)}}table of binned frequencies and variables used for plotting and diagnostics{p_end}
{synopt:{cmd:e(G)}}delta-method Jacobian for transformed parameters, when available{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Functions}{p_end}
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
