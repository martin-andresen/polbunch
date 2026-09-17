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
{synopt:{opt w:indow(numlist)}}two integers specifying the number of excluded bins below and above the cutoff; default is {cmd:window(1 0)}{p_end}
{synopt:{opt est:imator(#)}}estimator to use; default is {cmd:estimator(3)}{p_end}
{synopt:{opt delta:max(#)}}upper bound on the structural shift {cmd:delta} in the estimator 2/3 profile search; default {cmd:deltamax(1)} (a 100% earnings response){p_end}
{synopt:{opt log}}specifies that the running variable is in logs{p_end}
{synopt:{opt t0(#)}}linear tax rate below the cutoff{p_end}
{synopt:{opt t1(#)}}linear tax rate above the cutoff{p_end}

{syntab:Reporting and transformations}
{synopt:{opt notransform}}report raw estimating-equation coefficients rather than transformed bunching parameters{p_end}
{synopt:{opt constant}}use the constant-density approximation when turning bunching into a response and elasticity; default for estimators 1 and 2{p_end}
{synopt:{opt exact}}invert the counterfactual-density integral exactly when turning bunching into a response and elasticity; default for estimators 0 and 3{p_end}
{synopt:{opt poolmass}}form the bunching mass as M minus the integral of {cmd:h0} over the whole excluded region; default for estimators 1 and 2{p_end}
{synopt:{opt splitmass}}form the bunching mass as M minus the integral of {cmd:h0} up to the cutoff and of {cmd:h1} above it; default for estimators 0 and 3{p_end}

{syntab:Estimation controls}
{synopt:{opt allownegative}}let the structural shift {cmd:delta} range down to -1 in the estimator 2/3 profile; by default the search is restricted to {cmd:delta} >= 0{p_end}
{synopt:{opt positive}}restrict the structural shift to be positive (now the default for estimators 2/3; retained for backward compatibility){p_end}
{synopt:{opt nonormalize}}estimate on the original running-variable scale instead of normalizing bin centers by the cutoff and bin width{p_end}
{synopt:{opt nozero}}do not fill empty bins with zero counts; by default empty bins inside the observed support are included{p_end}
{synopt:{opt nodrop}}do not drop endpoint bins that appear to be cut by sample selection{p_end}
{synopt:{opt norankred}}do not reduce the polynomial degree when the unrestricted regression is rank-deficient; by default {cmd:polbunch} reduces the degree automatically. When {cmd:norankred} actually keeps a rank-deficient order the specification tests are disabled (they compare against the degenerate unrestricted fit){p_end}
{synopt:{opt norankcheck}}skip the separate-sides rank check entirely and estimate at the requested polynomial degree; also disables specification tests (the analytical bias is still reported){p_end}

{syntab:Inference}
{synopt:{opt vce(vcetype)}}{it:vcetype} is {cmd:conventional} (default; synonym {cmd:analytic}), {cmd:robust}, {cmd:hc2}, {cmd:hc3}, {cmd:cluster} {it:clustvar} (or {cmd:cluster} {it:stub}{cmd:, nclusters(}{it:#}{cmd:)} for pre-binned data), {cmd:bootstrap}{cmd:[, }{it:subopts}{cmd:]}, or {cmd:none}{p_end}
{synopt:{opt scale(spec)}}overdispersion multiplier for {cmd:vce(conventional)}, as in {helpb glm}; {it:spec} is {cmd:1} (default), {cmd:x2} (Pearson {it:phi}-hat), or a positive number{p_end}
{synopt:{opt reps(#)}}number of bootstrap repetitions for {cmd:vce(bootstrap)}; default {cmd:reps(500)}; also settable as {cmd:vce(bootstrap, reps(#))}{p_end}
{synopt:{opt nomasscorr}}disable the reference-region overdispersion correction to the bunching-mass row under {cmd:vce(robust)}/{cmd:hc2}/{cmd:hc3}{p_end}
{synopt:{opt nodots}}suppress bootstrap progress dots{p_end}
{synopt:{it:vce(bootstrap) subopts}}{cmd:multinomial} (default) {c |} {cmd:residual} {c |} {cmd:wild}; {cmd:bayesian}; {cmd:normal} (default) {c |} {cmd:bc} {c |} {cmd:percentile}; {cmd:wildweights(rademacher|mammen|webb)}; {cmd:reps(#)}; {cmd:seed(#)}{p_end}
{synopt:{opt savebins(file[, replace])}}with {cmd:vce(cluster} {it:clustvar}{cmd:)} on individual data: write the histogram + raw co-visitation matrix for later {cmd:vce(cluster} {it:stub}{cmd:)} use{p_end}

{syntab:Model-restriction test}
{synopt:{opt test(string)}}which model-restriction test to report: {cmd:all} (default for estimators 1--3), {cmd:minimumdistance} (omnibus, chi2 K+1), {cmd:hausman} (focused on the elasticity, chi2 1), {cmd:wald} (default for estimator 4), or {cmd:none}{p_end}
{synopt:{opt contrast}}paired-bootstrap test of this estimator's elasticity against the model-consistent efficient reference ({cmd:estimator(3)}, {cmd:exact}, {cmd:splitmass}), refit on every resample; requires {cmd:vce(bootstrap)}, {cmd:t0()}/{cmd:t1()} and {cmd:estimator(1/2/3)}. See {help polbunch##contrast:Elasticity contrast}{p_end}
{synoptline}
{p2colreset}{...}

{pstd}
{cmd:estat gof} after {cmd:polbunch} reports the reference-region goodness of fit -- deviance, Pearson {it:X-squared}, deviance {it:R-squared}, AIC/BIC and quasi-AIC. See {help polbunch##gof:Goodness of fit}.


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
The option {opt window(L H)} defines the excluded bunching region. If the cutoff lies inside a bin, that cutoff-crossing bin is included in the excluded region together with {it:L} bins below
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
for log earnings, it is additive in log earnings. Because {cmd:polbunch} handles only convex kinks, the marginal buncher reduces earnings toward the cutoff and the structural shift {cmd:delta} is non-negative by construction; the profile search enforces {cmd:delta} >= 0. Specify {opt allownegative} to open the search back to {cmd:delta} > -1 (useful only when the true shift is ~0 and finite-sample noise pulls the profile minimum slightly below zero).

{phang}
{cmd:estimator(4)} implements a Saez-style trapezoid approximation. It estimates a left reference height, a right reference height, and the bunching mass directly. This estimator does not use the nonlinear profile restrictions and does
not run the model-restriction test. With {opt notransform}, it reports {cmd:h0:_cons}, {cmd:h1:_cons}, and {cmd:bunching:B}. With the default transformation, it reports the reference heights, the number of bunchers, excess mass, and, when the trapezoid response equation can be solved, the shift, marginal response, and elasticity.

{pstd}
For estimators 2 and 3, {cmd:polbunch} uses a profile implementation of the corresponding nonlinear least-squares estimator. For any candidate value of the structural shift parameter, the remaining polynomial coefficients enter the model linearly
and are therefore estimated by least squares. The command then minimizes the resulting profiled sum of squared residuals over the structural parameter only. This recovers the same minimum and the same parameter estimates as the full nonlinear
least-squares problem, but avoids repeatedly optimizing over all polynomial coefficients jointly.
{p_end}

{pstd}
The one-dimensional profile minimization uses a coarse grid over {cmd:delta} followed by a golden-section refinement around every local minimum found, rather than a single Newton start. When the model restriction is misspecified
the profiled sum of squares can be multi-modal, so a single start can converge to different basins as the polynomial order or window changes; the grid search removes that dependence. {cmd:delta} is confined to {bf:[0, deltamax]} by default
({cmd:deltamax(1)} = a 100% earnings response at the kink) -- the non-negativity is implied by the convex kink, and it also removes a spurious degenerate branch near {cmd:delta} = -1. Raise {opt deltamax()} if a larger structural response is
plausible; specify {opt allownegative} to search {bf:(-1, deltamax]} instead.
{p_end}

{pstd}
If the profile minimizer lands on a bound of the search interval, or the profiled sum of squares has a second interior local minimum at a materially different {cmd:delta} whose fit is within 10% of the best, the structural shift is treated as weakly identified: {cmd:polbunch}
reports the counterfactual density, number of bunchers and excess mass (and, for estimator 2, the {cmd:delta} column), but withholds the shift, marginal response and elasticity, sets {cmd:e(hasresp)} to 0 and {cmd:e(delta_weakid)} to 1, and displays a
note. This most often means the excluded window is too narrow to pin down the response length or the level-shift restriction is rejected -- consult the minimum-distance or Hausman test and try a different {opt polynomial()} or {opt window()}.
{p_end}

{marker transform}{...}
{title:Transformed parameters}

{pstd}
Unless {opt notransform} is specified, {cmd:polbunch} transforms the raw estimating-equation parameters into economically interpretable quantities. The transformed output may include the estimated counterfactual density under the low-tax regime {cmd:h0}, the density under the high-tax regime {cmd:h1}, the number of bunchers, excess mass, the proportional shift, the response of the marginal buncher, and, if {opt t0()} and {opt t1()} are specified, the elasticity.

{pstd}
For estimator 3 under the default {opt splitmass} and {opt exact}, the response and elasticity are the model-implied structural shift {cmd:delta} itself. If {opt log} is specified, the response is a displacement in log earnings and the elasticity is the log response divided by the log net-of-tax-rate change. Without {opt log}, the level response is converted to a proportional shift before computing the elasticity. Under {opt poolmass}, estimator 3 instead backs the response out of the pooled reduced-form bunching mass, exactly as estimators 0, 1 and 2 do, and reports the structural {cmd:delta} in a separate column for reference.

{pstd}
For estimator 4, the Saez transformation first computes excess mass using the average of the two trapezoid endpoints. It then attempts to invert the trapezoid response equation. If the equation has no real positive solution, {cmd:polbunch} reports the reference heights, number of bunchers, and excess mass, but omits shift, marginal response, and elasticity and displays a note.

{pstd}
The {opt constant} / {opt exact} pair controls how bunching is converted to a response and elasticity. {opt constant} uses a constant-density approximation (the marginal response is the bunching mass divided by the counterfactual height at the cutoff); it can be useful for comparison with older procedures but may be biased when the density changes substantially over the response region. {opt exact} instead solves the counterfactual-density integral equation for the response. If neither is given, {cmd:polbunch} uses {opt constant} for estimators 1 and 2 and {opt exact} for estimators 0 and 3. The two are mutually exclusive.

{pstd}
The {opt poolmass} / {opt splitmass} pair controls how the bunching mass B is formed from the observed mass M in the excluded region. {opt poolmass} sets B = M minus the integral of the estimated counterfactual density {cmd:h0} over the entire excluded region. {opt splitmass} sets B = M minus the integral of {cmd:h0} from the lower limit to the cutoff, minus the integral of the post-tax density {cmd:h1} from the cutoff to the upper limit; this is the theoretically grounded calculation. If neither is given, {cmd:polbunch} uses {opt poolmass} for estimators 1 and 2 and {opt splitmass} for estimators 0 and 3. The choice has no effect for estimator 1 (where {cmd:h0} equals {cmd:h1}) and is ignored for estimator 4. The two are mutually exclusive.


{marker tests}{...}
{title:Model-restriction tests}

{pstd}
For estimators 1--4, {cmd:polbunch} tests the restrictions the selected estimator imposes on the earnings density across the kink, against the unrestricted two-sided fit ({cmd:estimator(0)} -- a separate degree-{it:K} polynomial on each side). Testing is skipped under {cmd:test(none)} or {cmd:vce(none)}. Two kinds of test are available, controlled by {opt test(string)}.

{phang}
{cmd:test(minimumdistance)} -- the {bf:omnibus} test. A minimum-distance / overidentification statistic on the whole cross-kink coefficient vector (the {it:K}+1 polynomial-coefficient jumps plus the excess mass), profiling out the estimator's structural nuisance parameter. Distributed chi-squared with {it:K}+1 degrees of freedom. For the naive estimator this coincides with {cmd:test(wald)} (there is no nuisance parameter), so it is not offered separately there. Available for estimators 2 and 3.

{phang}
{cmd:test(hausman)} -- the {bf:focused} test. A generalised (Wooldridge-form) Hausman statistic on the single parameter of interest: it contrasts the restricted elasticity with the elasticity implied by the unrestricted fit, standardised by the variance of the {it:difference} -- built from the joint influence functions of the two estimators on the shared bin counts, so their covariance is accounted for and no efficiency assumption on the restricted estimator is required. Distributed chi-squared with 1 degree of freedom. Requires an analytic {cmd:vce()} ({cmd:analytic}, {cmd:robust}, {cmd:hc2}, {cmd:hc3} or {cmd:cluster}); the meat matrix carries through. Available for estimators 1, 2 and 3. {cmd:polbunch} also reports the unrestricted elasticity and its standard error ({cmd:e(elast_unrestricted)}, {cmd:e(se_elast_unrestricted)}): a large standard error there -- wild extrapolation of the counterfactual into the excluded region -- means the test has little power in that application and the omnibus statistic is the one to read. {it:e_U} is always computed on the {cmd:exact} + {cmd:splitmass} axes -- the consistent combination -- regardless of the reported estimate's axes. So for the model-consistent estimator the contrast isolates the cross-kink density restriction; for the naive and Chetty estimators it additionally reflects their own {cmd:constant}/{cmd:poolmass} approximations, which is appropriate, since those are part of the estimator's bias. When the exact inversion has no real root on the extrapolated unrestricted counterfactual -- common with a high {opt polynomial()} and a wide excluded region -- {it:e_U} cannot be formed and the Hausman test is not reported; lower the degree, narrow {opt window()}, or use {cmd:minimumdistance}, which needs no inversion.

{phang}
{cmd:test(wald)} (default for estimators 1 and 4) is a Wald test of the restriction the selected estimator imposes, evaluated at the unrestricted estimates. For the naive estimator it is the omnibus test. When specified for estimators 2 or 3 it plugs in a mass-implied value of the structural parameter rather than profiling it, so {cmd:polbunch} prints a note pointing to {cmd:minimumdistance} as the formal test.

{phang}
{cmd:test(all)} runs every applicable test: {cmd:wald} + {cmd:hausman} for the naive estimator, {cmd:minimumdistance} + {cmd:hausman} for estimators 2 and 3, {cmd:wald} for Saez.

{phang}
{cmd:test(none)} suppresses all model-restriction testing.

{pstd}
{bf:What a rejection means.} Both tests check whether the density prediction of the estimator's structural model (assumption A1) holds, given the polynomial counterfactual (A3). They do {it:not} test the part of A1 both estimators share -- the marginal-buncher inversion that maps excess mass to an elasticity -- which cannot be tested from a single kink. A non-rejection therefore means "the assumed response is consistent with the shape of the density," not "the elasticity is correct." Where the omnibus and focused statistics agree, the detected misspecification loads on the elasticity; where they diverge (omnibus rejects, Hausman does not), the lack of fit is in a direction the elasticity is insensitive to.


{marker contrast}{...}
{title:Elasticity contrast}

{pstd}
{opt contrast} adds a paired-bootstrap comparison of the {it:reported} estimator's elasticity with the elasticity from the {bf:model-consistent efficient reference} -- {cmd:estimator(3)} with {cmd:exact} inversion and {cmd:splitmass} -- fitted on the {bf:same reference bins and window}. It is available for {cmd:estimator(1)}, {cmd:(2)} and {cmd:(3)} (with {cmd:constant} or {cmd:poolmass}), requires {cmd:vce(bootstrap[, ...])} and {opt t0()}/{opt t1()}, and errors when the reported estimate already {it:is} the reference.

{pstd}
On every bootstrap replication {cmd:polbunch} refits the reference on the {it:same} resampled histogram the reported estimator saw, so the two elasticity draws are paired. The bootstrap SD of the difference {it:d* = e_reported* - e_reference*} is its standard error {bf:with the two estimators' dependence built in} -- there is no efficiency assumption and no analytic covariance to derive, and it works for every {cmd:vce(bootstrap)} flavour ({cmd:multinomial}, {cmd:residual}, {cmd:wild}, {cmd:bayesian}). A replication on which the reference profile is weakly identified or the exact inversion has no real root contributes a missing value and is dropped pairwise; the count of usable pairs is reported.

{pstd}
The block reports both point elasticities, their difference, the paired bootstrap SE, {cmd:corr(reported, reference)} across replications (typically well above 0.9 -- the two estimators share the polynomial fit -- which is precisely why quadrature of the two marginal SEs is far too conservative), a 95% normal-approximation CI, a {it:z} statistic with its normal p-value, and a bootstrap p-value (the share of replications with {bf:|d* - d| >= |d|}).

{pstd}
{bf:Interpretation.} The reference is consistent for the structural elasticity under A1 and the common counterfactual; the reported estimator generally is not (its constant-density inversion, mass pooling, or -- for {cmd:estimator(1)} -- no counterfactual adjustment at all, each carry a probability limit of their own). So {it:d} measures how far this estimator's approximation and finite-sample behaviour move the answer {it:on this sample}, relative to sampling noise -- it is not a signed statement that one estimator is biased, and it is {bf:not comparable across estimators}. Because {it:corr} is large the SE of {it:d} is small and the test has high power: it will flag even economically minor gaps (for {cmd:estimator(1)}/{cmd:(2)} against {cmd:estimator(3)}, almost always). Read the reported difference for magnitude alongside the p-value for significance. When it does {it:not} reject -- e.g. {cmd:estimator(3) constant} vs {cmd:exact} on a near-symmetric window with a small response -- that is genuine reassurance that the approximation is immaterial in that application.

{pstd}
Results are stored in {cmd:e(contrast_diff)}, {cmd:e(contrast_se)}, {cmd:e(contrast_corr)}, {cmd:e(contrast_z)}, {cmd:e(contrast_p)}, {cmd:e(contrast_p_pctile)}, {cmd:e(contrast_ci_ll)}, {cmd:e(contrast_ci_ul)}, {cmd:e(contrast_elast)}, {cmd:e(contrast_elast_ref)}, {cmd:e(contrast_reps)}, {cmd:e(contrast_reps_used)} and {cmd:e(contrast_ref)}. After the command {cmd:e(b)} is unchanged, so {cmd:lincom} and {cmd:test} still refer to the reported model.

{pstd}
For a contrast between {it:any} two {cmd:polbunch} specifications -- a different estimator, polynomial order, inversion or mass axis, or window, and estimands other than the elasticity -- fit and {helpb estimates store} both, then use the postestimation command {helpb polbunch_contrast}, which resamples the shared histogram (read from {cmd:e(bins)}) and re-fits both specifications on every replication.


{marker bias}{...}
{title:Analytical bias}

{pstd}
Unless {opt nobias} is specified, and when {opt t0()} and {opt t1()} are given, {cmd:polbunch} appends the analytical bias of the fitted estimator (computed by {helpb polbunchbias} in post-estimation mode) to the transformed-parameter table and stores it in {cmd:e(bias_elasticity)}, {cmd:e(bias_B)}, {cmd:e(bias_shift)} and related scalars.

{pstd}
This bias is measured against {ul:the counterfactual that the selected estimator itself assumes} -- the fitted degree-K polynomial for estimators 0--3, the Saez two-point line for estimator 4. It is an internal-consistency diagnostic, {bf:not} a quantity that can be compared across the {opt estimator()} settings: each estimator posits a different counterfactual, so a small bias figure for one and a large one for another does not rank them. Estimators 0 and 3 (under {opt exact} + {opt splitmass}) and estimator 4 (under {opt constant}) are internally consistent by construction, so no bias line is shown for them. For a comparison across estimators, run {helpb polbunchbias} standalone with the same {opt h0poly()} counterfactual supplied to each.


{marker gof}{...}
{title:Goodness of fit}

{pstd}
For every analytic {cmd:vce()} (and {cmd:vce(bootstrap, residual|wild)}), {cmd:polbunch} evaluates how well the fitted counterfactual describes the {it:reference} bins -- {bf:h0} below the cutoff and, for the restricted estimators, {bf:h1(delta-hat)} above it. The excess-mass bins are not part of this fit and do not enter. The histogram is treated as multinomial in {it:N} and, conditional on {it:N}, independent Poisson, so the diagnostics are the usual Poisson / quasi-Poisson ones:

{p 8 12 2}{cmd:e(deviance)}, {cmd:e(pearson_x2)} -- absolute fit; each {it:approx chi-squared(e(gof_df))} under correct specification. {cmd:e(pearson_x2)}/{cmd:e(gof_df)} is exactly {cmd:e(dispersion)}. Their p-values ({cmd:e(deviance_p)}, {cmd:e(pearson_x2_p)}) reject for essentially any real bunching dataset because {it:N} is huge; read the {it:sizes}, not the tests.{p_end}
{p 8 12 2}{cmd:e(r2_dev)} -- deviance R-squared (Cameron and Windmeijer 1997), {bf:1 - deviance/deviance_null}, in [0,1]. The share of the null (intercept-only) deviance explained by the counterfactual.{p_end}
{p 8 12 2}{cmd:e(aic)}, {cmd:e(bic)}, {cmd:e(qaic)}, {cmd:e(qaicc)} -- information criteria for choosing among {opt estimator()} settings and {opt polynomial()} degrees. Use the {it:quasi} versions {cmd:e(qaic)}/{cmd:e(qaicc)} ({bf:-2*loglik/phi-hat + 2p}) when {cmd:e(dispersion)} > 1, which it usually is. For a formal ranking across a set of candidate models, recompute QAIC by hand with one common {it:phi}-hat -- conventionally the most general model in the set, e.g. {cmd:estimator(0)} at the highest degree ({cmd:e(gof_ll)}, {cmd:e(gof_np)} and {cmd:e(gof_nbins)} are stored for this).{p_end}

{pstd}
The fit is also reported {it:split at the cutoff} ({cmd:e(dispersion_below)}, {cmd:e(dispersion_above)}). Below {it:z*} nobody bunches and {it:h1 = h0}, so overdispersion there is a near-pure test of A3 (the polynomial's fit to {it:h0}); above {it:z*} it also picks up an A1/A2 shape error in the response-implied {it:h1}. Combined with a sweep of {opt polynomial()}: if {cmd:e(dispersion_below)} falls toward 1 as the degree rises, the misfit was polynomial-approximation error (A3, curable); if it stays above 1, the counterfactual has non-polynomial structure (a second kink, heaping) that more degrees will not fix; {cmd:e(dispersion_below)} near 1 with the restriction test still rejecting points instead to the behavioural model (A1). Under panel or repeated observations run a single cross-section or {cmd:vce(cluster)} first, so a raised dispersion is misfit and not clustering.

{pstd}
A compact summary line is printed under the coefficient table; {cmd:estat gof} shows the full block and returns it in {cmd:r()}. For {opt estimator(4)} the counterfactual is deliberately a two-level local approximation, so its reference-region fit will look poor and should not be compared head-to-head with the polynomial estimators.


{marker inference}{...}
{title:Inference}

{pstd}
The variance estimator is controlled by {opt vce(vcetype)}, following the usual Stata form {cmd:vce(}{it:type}{cmd:[, }{it:subopts}{cmd:])}. All types leave the point estimates (and the specification tests' point estimates) unchanged; they differ only in the standard errors, the covariance matrix and the test statistics' weight matrices.

{phang}
{cmd:vce(conventional)} (the default; synonyms {cmd:analytic}, {cmd:unadjusted}, {cmd:oim}) computes the collapsed-data delta-method standard errors from the model-implied multinomial variance of the bin counts, {bf:Var(y_j) ~ y_j}. It reproduces the variance of a regression run in the expanded one-row-per-individual-per-bin data, without building that data, and is the same object as Saez (2010)'s delta method. Efficient and correct when the counts really are multinomial (independent individuals, correct polynomial).

{phang}
{cmd:scale(x2)} multiplies the {cmd:vce(conventional)} covariance by the Pearson overdispersion {bf:phi-hat = sum_j (y_j - yhat_j)^2 / yhat_j / (n_ref - p)} estimated from the reference bins -- the quasi-Poisson standard error, as in {helpb glm}. {cmd:scale(#)} uses a literal factor. Supported only with {cmd:vce(conventional)}; {cmd:phi-hat} is returned in {cmd:e(dispersion)} for every analytic run regardless.

{phang}
{cmd:vce(robust)} (= {cmd:hc1}; {cmd:vce(hc0)}/{cmd:vce(hc2)}/{cmd:vce(hc3)} select the other finite-sample corrections) replaces the multinomial diagonal with the Eicker-White residual meat {bf:(y_j - yhat_j)^2}, HC-corrected. The bunching-mass row, whose fitted residual is ~0 by construction, is instead lifted to {bf:phi-hat * Hstar(1 - Hstar/N)} -- borrowing the reference-region overdispersion into the excluded window (turn this off with {opt nomasscorr}). Robust to arbitrary bin-level heteroskedasticity, to polynomial misspecification treated as noise (round-number heaping, secondary bumps, a neighbouring kink), and to the per-bin part of overdispersion from repeated individuals or year effects. {it:Not} robust to cross-bin correlation (the off-diagonals stay multinomial). Typically wider than {cmd:vce(conventional)} and close to the residual-bootstrap standard errors common in the bunching literature. Replacing only the diagonal of an otherwise-multinomial meat is not guaranteed to leave it positive semi-definite; when it isn't, the reported standard errors come from the nearest valid covariance instead (see {cmd:e(vce_psdclip)}).

{phang}
{cmd:vce(bootstrap)} re-runs the full estimator (integration constraint and all) on {opt reps(#)} resampled bin-count vectors and reports the standard deviation of the resulting estimates. The {it:subopts} select the resampling scheme:

{phang2}{cmd:multinomial} (default) -- Dirichlet resample of the bin counts. Equivalent to resampling individuals and re-binning, so it targets the same object as {cmd:vce(conventional)}.{p_end}
{phang2}{cmd:residual} -- resample, {it:iid} with replacement, the reference-bin residuals of the {it:fitted model itself} -- {bf:y_j - h0(z_j)} below the kink and {bf:y_j - h1(z_j; delta-hat)} above -- and add them back to that same fitted mean (Chetty/CFOP). Because the pool is drawn around the estimated counterfactual (not a separate, more flexible per-side fit), it carries the model's own lack of fit; it imposes a common residual variance across bins, so it is close to {cmd:scale(x2)} in spirit. The integration constraint is re-solved each draw. The excluded bins carry no per-bin counterfactual, so they are perturbed by {bf:sqrt(phi-hat * y_j)} noise -- the resampling twin of the {cmd:vce(robust)} mass row.{p_end}
{phang2}{cmd:wild} -- multiply each reference-bin residual (of the fitted model, as above) by a mean-0 variance-1 weight ({cmd:wildweights(rademacher|mammen|webb)}, default {cmd:rademacher}). Heteroskedasticity-robust: the resampling twin of {cmd:vce(robust)}, and the appropriate version of the residual bootstrap when densities heap.{p_end}
{phang2}{cmd:bayesian} -- Bayesian bootstrap (Dirichlet weights); like {cmd:multinomial}, targets the individual-sampling variance.{p_end}

{pstd}
The displayed coefficient table always shows normal-approximation confidence intervals. {cmd:bc} or {cmd:percentile} additionally store bias-corrected or percentile bootstrap intervals in {cmd:e(ci_bc)} / {cmd:e(ci_percentile)}. {cmd:reps(#)} and {cmd:seed(#)} may also be given inside {cmd:vce(bootstrap, ...)}.

{phang}
{cmd:vce(cluster} {it:clustvar}{cmd:)} computes cluster-robust standard errors for pooled panel data, where the same individuals recur across years. Before binning, {cmd:polbunch} forms the bin co-visitation matrix {bf:M = sum_i c_i c_i'}, where {bf:c_i} is cluster {it:i}'s vector of visit counts over the rows of the stacked estimating system, and uses the cluster-robust meat {bf:n/(n-1) (M - y y'/n)} in the sandwich. It reduces to {cmd:vce(conventional)} (up to {bf:n/(n-1)}) when no individual appears more than once, and is the analytic counterpart of a bootstrap that resamples individuals -- far cheaper on large panels. It picks up both the within-individual over-dispersion and the positive cross-bin correlation that {cmd:vce(robust)} misses. The number of clusters is returned in {cmd:e(N_clust)}.

{phang2}
With {it:individual-level} input, {it:clustvar} is a cluster identifier and {bf:M} is built internally.{p_end}

{phang2}
With {it:pre-binned} input, {it:clustvar} is instead a {it:stub} naming {bf:J} variables {it:stub}{cmd:1}..{it:stub}{cmd:J} that hold the {bf:J x J} bin co-visitation matrix, one column per bin, row-aligned with the histogram; the number of clusters must be given as {cmd:vce(cluster} {it:stub}{cmd:, nclusters(}{it:#}{cmd:))} (or carried in the dataset characteristic {cmd:_dta[polbunch_nclusters]}). The matrix must be symmetric with {bf:M[j,j] = sum_i c_ij^2 >= y_j}. Only the full matrix is accepted -- a diagonal-only {bf:M} does not give a valid covariance; use {cmd:vce(robust)} on the histogram as the per-bin-only hedge.{p_end}

{phang}
{opt savebins(filename[, replace])} (individual data, with {cmd:vce(cluster} {it:clustvar}{cmd:)}) writes the histogram and the raw {bf:J x J} co-visitation matrix to {it:filename}{cmd:.dta} as {cmd:freq}, {cmd:midpoint}, {cmd:cv1}..{cmd:cvJ}, with the cluster count, cutoff and bandwidth stored as dataset characteristics. The file is exactly the input {cmd:vce(cluster cv)} expects, so a panel can be reduced to a disclosable histogram-plus-matrix once and analysed later without the microdata.

{phang}
{cmd:vce(none)} suppresses all internal variance estimation -- point estimates only. Useful in Monte Carlo work or under Stata's {cmd:bootstrap} prefix.

{pstd}
{cmd:vce(conventional)}, {cmd:vce(robust)} and the bootstrap types assume the bin counts are independent across bins and the bandwidth and binning scheme are fixed. Clustering from repeated observations of the same individual across pooled years is handled by {cmd:vce(cluster} {it:clustvar}{cmd:)} -- from the microdata, or from a pre-binned histogram plus the supplied co-visitation matrix; other higher-level clustering requires a resampling procedure outside {cmd:polbunch}.


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
{opt norankred} suppresses the automatic polynomial degree reduction that {cmd:polbunch} performs when the unrestricted two-sided polynomial regression is rank-deficient. By default, the degree is reduced one step at a time until identification is restored; a note is displayed. The unrestricted fit is still run once (to seed starting values) and the specification tests, which compare the restricted estimator against that unrestricted fit, are still reported.

{phang}
{opt norankcheck} goes further: it skips the separate-sides identification check altogether and estimates with exactly the requested {opt polynomial()} degree. This is useful because a restricted estimator (2 or 3) imposes enough structure across the cutoff to be identified at a degree where two free one-sided polynomials are not. Because the specification tests compare the restricted estimator against an identified unrestricted fit, {opt norankcheck} disables them (equivalent to adding {cmd:test(none)}); a note is displayed. The analytical bias is unaffected -- it is read from the fitted counterfactual polynomial, tax rates, window and elasticity, not from the separate one-sided polynomials -- and is still reported. Recommended only with estimators 0, 2 and 3 -- estimator 1 is itself the unrestricted two-sided fit and gains nothing from forcing the order.

{phang}
{opt scale(spec)} multiplies the {cmd:vce(conventional)} covariance by an overdispersion factor, exactly as in {helpb glm}. {cmd:scale(x2)} uses the Pearson dispersion estimated from the reference bins; {cmd:scale(#)} uses a literal factor; {cmd:scale(1)} (the default) leaves the covariance unchanged. It is an error to combine {cmd:scale()} with anything other than {cmd:vce(conventional)}. The Pearson dispersion is always returned in {cmd:e(dispersion)} for analytic runs whether or not {cmd:scale()} is used.

{phang}
{opt nomasscorr} turns off the reference-region overdispersion correction that {cmd:vce(robust)}/{cmd:hc2}/{cmd:hc3} apply to the bunching-mass row. With the correction on (the default), that row's variance is {bf:phi-hat * Hstar(1 - Hstar/N)} rather than the pure-Poisson {bf:Hstar(1 - Hstar/N)}; the former makes the analytic robust standard error agree with the wild bootstrap and matters most for estimator 1, where the excess mass is dominated by the observed mass rather than by the extrapolated counterfactual.

{phang}
{opt test(string)} selects the restriction test. The default is {cmd:all} for estimators 2 and 3 and {cmd:wald} for estimators 1 and 4. Specifying {cmd:test(none)} skips all testing. Note that testing is also suppressed when {cmd:vce(none)} is in effect.


{marker examples}{...}
{title:Examples}

{pstd}
Generate simulated data from the companion data-generating command:{p_end}

{phang2}{cmd:. polbunchgendata, obs(10000) t0(0.2) t1(0.6) elasticity(0.4) cutoff(1)}{p_end}

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
Use the bootstrap for inference with 200 repetitions:{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(1) vce(bootstrap) reps(200)}{p_end}

{pstd}
Report misspecification-robust (Eicker-White residual) standard errors, comparable to the residual-bootstrap standard errors common in the bunching literature:{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(7) t0(0.2) t1(0.6) vce(robust)}{p_end}

{pstd}
Quasi-Poisson standard errors (one overdispersion parameter):{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(7) t0(0.2) t1(0.6) scale(x2)}{p_end}

{pstd}
Reproduce the Chetty/CFOP residual bootstrap, or its heteroskedasticity-robust (wild) counterpart with percentile intervals:{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(7) t0(0.2) t1(0.6) vce(bootstrap, residual)}{p_end}
{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(7) t0(0.2) t1(0.6) vce(bootstrap, wild percentile) reps(999)}{p_end}

{pstd}
Test whether the Chetty estimate's elasticity differs significantly from the model-consistent efficient estimator, accounting for the two being estimated on the same data:{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(7) t0(0.2) t1(0.6) estimator(2) vce(bootstrap) reps(999) contrast}{p_end}

{pstd}
Suppress internal variance estimation, for example when using Stata's bootstrap prefix:{p_end}

{phang2}{cmd:. bootstrap, reps(200): polbunch z, cutoff(1) bw(0.01) polynomial(1) vce(none)}{p_end}

{pstd}
Report only the focused Hausman test on the elasticity:{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(1) test(hausman)}{p_end}

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
{synopt:{cmd:e(constant)}}1 if the constant-density response inversion was used, 0 if the exact inversion was used{p_end}
{synopt:{cmd:e(nosplit)}}1 if the bunching mass used the {opt poolmass} calculation, 0 if it used {opt splitmass}{p_end}
{synopt:{cmd:e(chi2_wald)}}Wald test chi-squared statistic, when computed{p_end}
{synopt:{cmd:e(p_wald)}}Wald test p-value, when computed{p_end}
{synopt:{cmd:e(df_wald)}}Wald test degrees of freedom, when computed{p_end}
{synopt:{cmd:e(chi2_minimumdistance)}}minimum-distance test chi-squared statistic, when computed{p_end}
{synopt:{cmd:e(p_minimumdistance)}}minimum-distance test p-value, when computed{p_end}
{synopt:{cmd:e(df_minimumdistance)}}minimum-distance test degrees of freedom, when computed{p_end}
{synopt:{cmd:e(chi2_hausman)}}focused Hausman test chi-squared statistic (elasticity contrast, 1 df), when computed{p_end}
{synopt:{cmd:e(p_hausman)}}focused Hausman test p-value, when computed{p_end}
{synopt:{cmd:e(df_hausman)}}focused Hausman test degrees of freedom (1), when computed{p_end}
{synopt:{cmd:e(elast_unrestricted)}, {cmd:e(se_elast_unrestricted)}}elasticity implied by the unrestricted two-sided fit and its standard error (reported with the Hausman test; a large standard error flags weak power){p_end}
{synopt:{cmd:e(delta_md)}}structural shift estimate from the minimum-distance test, when computed{p_end}
{synopt:{cmd:e(dispersion)}}Pearson overdispersion {it:phi}-hat from the reference bins, for analytic {cmd:vce()}{p_end}
{synopt:{cmd:e(deviance)}}Poisson deviance of the counterfactual on the reference bins ({cmd:e(deviance_p)} its {it:chi-squared} p-value); {cmd:e(deviance_null)} for the intercept-only fit. See {help polbunch##gof:Goodness of fit}{p_end}
{synopt:{cmd:e(pearson_x2)}}Pearson {it:X-squared} on the reference bins ({cmd:e(pearson_x2_p)} its p-value); {cmd:e(pearson_x2)} / {cmd:e(gof_df)} equals {cmd:e(dispersion)}{p_end}
{synopt:{cmd:e(r2_dev)}}deviance {it:R-squared} (Cameron-Windmeijer), {bf:1 - e(deviance)/e(deviance_null)}{p_end}
{synopt:{cmd:e(gof_ll)}}Poisson log-likelihood of the counterfactual on the reference bins ({cmd:e(gof_ll_null)} for the intercept-only fit){p_end}
{synopt:{cmd:e(aic)}, {cmd:e(bic)}}Akaike / Bayesian information criteria from {cmd:e(gof_ll)} and {cmd:e(gof_np)}{p_end}
{synopt:{cmd:e(qaic)}, {cmd:e(qaicc)}}quasi-AIC and its small-sample-corrected form, {bf:-2*ll/phi-hat + 2p} -- for comparing estimators / degrees under overdispersion{p_end}
{synopt:{cmd:e(gof_nbins)}, {cmd:e(gof_np)}, {cmd:e(gof_df)}}reference bins entering the fit, fitted parameters, and their difference{p_end}
{synopt:{cmd:e(gof_rmse)}}root mean squared reference-bin residual{p_end}
{synopt:{cmd:e(gof_massresid)}}relative residual of the mass-restriction row ({it:approx} 0; a large value for estimator 2/3 flags an unconverged {it:delta} solve){p_end}
{synopt:{cmd:e(dispersion_below)}, {cmd:e(dispersion_above)}}Pearson dispersion of the reference fit {it:below} and {it:above} the cutoff (with {cmd:e(deviance_below)} / {cmd:e(deviance_above)} and their df). Below the cutoff there is no bunching and {it:h1 = h0}, so {cmd:e(dispersion_below)} is a near-pure test of A3; {cmd:e(dispersion_above)} also reflects an A1/A2 shape error in {it:h1}. See {help polbunch##gof:Goodness of fit}{p_end}
{synopt:{cmd:e(scalefactor)}}factor actually applied by {opt scale()}, when not 1{p_end}
{synopt:{cmd:e(masscorr)}}1 if the bunching-mass overdispersion correction was applied ({cmd:vce(robust)}/{cmd:hc2}/{cmd:hc3}){p_end}
{synopt:{cmd:e(vce_psdclip)}}{cmd:vce(robust)}/{cmd:hc2}/{cmd:hc3} only: 1 if the sandwich was not positive semi-definite and was projected onto the nearest valid covariance before standard errors were computed, 0 otherwise. Swapping only the diagonal of the multinomial meat for squared residuals is not guaranteed to leave it PSD; without this projection Stata's own {cmd:ereturn post} would silently report an all-zero variance (with only a terse, easily-missed warning) rather than error. A 1 here does not mean the fit is wrong -- it means the reported standard errors were the nearest valid ones to an otherwise-indefinite sandwich, most often when the polynomial order was reduced for multicollinearity or the sample sits far from the reference bins' center.{p_end}
{synopt:{cmd:e(bootreps)}}number of bootstrap repetitions, for {cmd:vce(bootstrap)}{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters, for {cmd:vce(cluster)}{p_end}
{synopt:{cmd:e(hasresp)}}1 if the shift, marginal response and elasticity are identified and reported; 0 if they were withheld{p_end}
{synopt:{cmd:e(delta_weakid)}}estimators 2/3: 1 if the structural shift {cmd:delta} was weakly identified (search bound hit, or multi-modal profile); 0 otherwise{p_end}
{synopt:{cmd:e(delta_nbasin)}}estimators 2/3: number of competing interior basins in the profile (1 = well identified; >= 2 triggers {cmd:e(delta_weakid)}){p_end}
{synopt:{cmd:e(delta_nonneg)}}estimators 2/3: 1 if the profile search was restricted to {cmd:delta} >= 0 (the default), 0 if {opt allownegative} was specified{p_end}
{synopt:{cmd:e(bias_elasticity)}}analytical bias of the fitted elasticity under the isoelastic model, {ul:measured against this estimator's own counterfactual and not comparable across {opt estimator()} settings} (and {cmd:e(bias_shift)}, {cmd:e(bias_B)}, {cmd:e(bias_slope)}, {cmd:e(bias_lambda)}, ...); see {help polbunch##bias:Analytical bias} and {help polbunchbias}{p_end}
{synopt:{cmd:e(bias_polynomial)}}polynomial order at which the analytical bias was evaluated{p_end}
{synopt:{cmd:e(bias_polyfull)}}requested polynomial order; when it exceeds {cmd:e(bias_polynomial)} the fitting design was too ill-conditioned for a bias at the full order and {cmd:polbunchbias} stepped down (a note is shown){p_end}
{synopt:{cmd:e(bias_l2err)}}relative change in the counterfactual from that order fallback; 0 when none{p_end}
{synopt:{cmd:e(bias_uniter_polynomial)}}set when the analytical-bias {cmd:iterate} loop failed to converge at this order and {cmd:e(bias_polynomial)} was PROMOTED (the default) to a lower order where it does converge: the demoted, un-iterated (plug-in) order{p_end}
{synopt:{cmd:e(bias_uniter_elasticity)}}(and {cmd:e(bias_uniter_shift)}, {cmd:e(bias_uniter_B)}, {cmd:e(bias_uniter_response)}) the demoted plug-in bias, kept alongside -- not discarded by -- the promoted {cmd:e(bias_elasticity)}; see {help polbunchbias}{p_end}
{synopt:{cmd:e(bias_iter_polynomial)}}with {cmd:nopromote} on the underlying {cmd:polbunchbias} call: the converged lower order, reported alongside (not replacing) {cmd:e(bias_elasticity)} instead of promoted into it{p_end}
{synopt:{cmd:e(bias_iter_elasticity)}}(and {cmd:e(bias_iter_shift)}, {cmd:e(bias_iter_B)}, {cmd:e(bias_iter_response)}) the {cmd:nopromote} converged-lower-order bias{p_end}
{synopt:{cmd:e(contrast_diff)}}with {opt contrast}: reported elasticity minus the reference ({cmd:estimator(3)}, {cmd:exact}, {cmd:splitmass}) elasticity{p_end}
{synopt:{cmd:e(contrast_se)}}paired-bootstrap standard error of {cmd:e(contrast_diff)} (carries the two estimators' dependence){p_end}
{synopt:{cmd:e(contrast_corr)}}bootstrap correlation between the reported and reference elasticity draws{p_end}
{synopt:{cmd:e(contrast_z)}, {cmd:e(contrast_p)}}{it:z} = {cmd:e(contrast_diff)}/{cmd:e(contrast_se)} and its two-sided normal p-value{p_end}
{synopt:{cmd:e(contrast_p_pctile)}}bootstrap p-value: share of replications with {bf:|d* - d| >= |d|}{p_end}
{synopt:{cmd:e(contrast_ci_ll)}, {cmd:e(contrast_ci_ul)}}95% normal-approximation CI for the difference{p_end}
{synopt:{cmd:e(contrast_elast)}, {cmd:e(contrast_elast_ref)}}the reported and reference elasticity point estimates{p_end}
{synopt:{cmd:e(contrast_reps)}, {cmd:e(contrast_reps_used)}}bootstrap replications requested, and those yielding a paired difference{p_end}
{p2colreset}{...}

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
{synopt:{cmd:e(vce)}}{cmd:vce()} type: {cmd:conventional}, {cmd:robust}/{cmd:hc0}..{cmd:hc3}, {cmd:cluster}, {cmd:bootstrap}, or {cmd:none}{p_end}
{synopt:{cmd:e(clustvar)}}cluster variable, for {cmd:vce(cluster} {it:clustvar}{cmd:)} on individual data{p_end}
{synopt:{cmd:e(covisstub)}}co-visitation stub, for {cmd:vce(cluster} {it:stub}{cmd:)} on pre-binned data{p_end}
{synopt:{cmd:e(boottype)}}bootstrap scheme: {cmd:multinomial}, {cmd:residual}, {cmd:wild}, or {cmd:bayesian}{p_end}
{synopt:{cmd:e(contrast_ref)}}with {opt contrast}: the reference specification, {cmd:estimator(3) exact splitmass}{p_end}
{synopt:{cmd:e(bootci)}}bootstrap CI type: {cmd:normal}, {cmd:bc}, or {cmd:percentile}{p_end}
{synopt:{cmd:e(wildweights)}}wild-bootstrap weight distribution, for {cmd:vce(bootstrap, wild)}{p_end}
{synopt:{cmd:e(scale)}}{opt scale()} specification, when not 1{p_end}
{synopt:{cmd:e(properties)}}usually {cmd:b V} when a variance matrix is posted{p_end}
{p2colreset}{...}

{synoptset 26 tabbed}{...}
{p2col 5 26 30 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector. With default transformation this contains density and bunching parameters; with {opt notransform} it contains raw estimating-equation coefficients{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of {cmd:e(b)}, when computed{p_end}
{synopt:{cmd:e(table)}}table of binned frequencies and variables used for plotting and diagnostics{p_end}
{synopt:{cmd:e(bins)}}the raw histogram actually fitted: bin count ({cmd:freq}) and bin {cmd:midpoint} in the original {it:z} units.  Read by {helpb polbunch_contrast} to refit the specification on resampled histograms{p_end}
{synopt:{cmd:e(G)}}delta-method Jacobian for transformed parameters, when available{p_end}
{synopt:{cmd:e(ci_bc)} / {cmd:e(ci_percentile)}}bias-corrected / percentile bootstrap confidence bounds for {cmd:e(b)}, when {cmd:vce(bootstrap, bc)} or {cmd:vce(bootstrap, percentile)} is used{p_end}
{p2colreset}{...}

{synoptset 26 tabbed}{...}
{p2col 5 26 30 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks the estimation sample{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{phang}
Andresen, Martin E. (2026). "A better polynomial bunching estimator", working paper.

{phang}
Cameron, A. Colin, and Frank A. G. Windmeijer (1997). "An R-squared measure of goodness of fit for some common nonlinear regression models", {it:Journal of Econometrics}.

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
