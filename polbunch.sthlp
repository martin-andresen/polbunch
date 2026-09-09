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
{synopt:{opt boot:reps(#)}}number of bootstrap repetitions for {cmd:vce(bootstrap)}; default {cmd:bootreps(500)}; also settable as {cmd:vce(bootstrap, reps(#))}{p_end}
{synopt:{opt nomasscorr}}disable the reference-region overdispersion correction to the bunching-mass row under {cmd:vce(robust)}/{cmd:hc2}/{cmd:hc3}{p_end}
{synopt:{opt nodots}}suppress bootstrap progress dots{p_end}
{synopt:{it:vce(bootstrap) subopts}}{cmd:multinomial} (default) {c |} {cmd:residual} {c |} {cmd:wild}; {cmd:bayesian}; {cmd:normal} (default) {c |} {cmd:bc} {c |} {cmd:percentile}; {cmd:wildweights(rademacher|mammen|webb)}; {cmd:reps(#)}; {cmd:seed(#)}{p_end}
{synopt:{opt savebins(file[, replace])}}with {cmd:vce(cluster} {it:clustvar}{cmd:)} on individual data: write the histogram + raw co-visitation matrix for later {cmd:vce(cluster} {it:stub}{cmd:)} use{p_end}

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
note. This most often means the excluded window is too narrow to pin down the response length or the level-shift restriction is rejected -- consult the minimum-distance or Hausman test and try a different {opt polynomial()} or {opt limits()}.
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
For estimators 1–4, {cmd:polbunch} tests the restrictions implied by the selected estimator, unless {cmd:test(none)} is specified or {cmd:vce(none)} is used. The test type is controlled by {opt test(string)}.

{phang}
{cmd:test(all)} (default for estimators 2 and 3) runs all applicable tests for the selected estimator and reports those that succeed. For estimators 2 and 3 this means minimum-distance and Hausman; for estimators 1 and 4 it reduces to the Wald test only.

{phang}
{cmd:test(hausman)} is a Hausman-type test that compares the restricted and unrestricted estimates. It uses the analytic variance of the difference between the two estimates when {cmd:vce(analytic)}, {cmd:vce(robust)} or {cmd:vce(cluster)} is in effect (the corresponding meat matrix carries through to the test's weight matrix), and the bootstrap covariance of the difference otherwise. Available for estimators 2 and 3.

{phang}
{cmd:test(wald)} (default for estimators 1 and 4) is a Wald test of the linear or nonlinear restrictions imposed by the selected estimator against the unrestricted estimator-0 estimates. For estimators 1 and 4 this reduces to a standard linear restriction test. When specified explicitly for estimators 2 or 3, {cmd:polbunch} prints a note recommending minimum-distance or Hausman instead, as the Wald statistic is a conditional shape diagnostic rather than a formal overall specification test.

{phang}
{cmd:test(minimumdistance)} is a minimum-distance test. It minimizes the Wald criterion over the structural parameter delta and compares the resulting minimum distance statistic to a chi-squared distribution. Available for estimators 2 and 3.

{phang}
{cmd:test(none)} suppresses all model-restriction testing.

{pstd}
All tests are reported as chi-squared statistics with associated p-values. They should be interpreted jointly as a test of the estimator's structural restrictions and the polynomial approximation used for the counterfactual density.


{marker bias}{...}
{title:Analytical bias}

{pstd}
Unless {opt nobias} is specified, and when {opt t0()} and {opt t1()} are given, {cmd:polbunch} appends the analytical bias of the fitted estimator (computed by {helpb polbunchbias} in post-estimation mode) to the transformed-parameter table and stores it in {cmd:e(bias_elasticity)}, {cmd:e(bias_B)}, {cmd:e(bias_shift)} and related scalars.

{pstd}
This bias is measured against {ul:the counterfactual that the selected estimator itself assumes} -- the fitted degree-K polynomial for estimators 0--3, the Saez two-point line for estimator 4. It is an internal-consistency diagnostic, {bf:not} a quantity that can be compared across the {opt estimator()} settings: each estimator posits a different counterfactual, so a small bias figure for one and a large one for another does not rank them. Estimators 0 and 3 (under {opt exact} + {opt splitmass}) and estimator 4 (under {opt constant}) are internally consistent by construction, so no bias line is shown for them. For a comparison across estimators, run {helpb polbunchbias} standalone with the same {opt h0poly()} counterfactual supplied to each.


{marker inference}{...}
{title:Inference}

{pstd}
The variance estimator is controlled by {opt vce(vcetype)}, following the usual Stata form {cmd:vce(}{it:type}{cmd:[, }{it:subopts}{cmd:])}. All types leave the point estimates (and the specification tests' point estimates) unchanged; they differ only in the standard errors, the covariance matrix and the test statistics' weight matrices.

{phang}
{cmd:vce(conventional)} (the default; synonyms {cmd:analytic}, {cmd:unadjusted}, {cmd:oim}) computes the collapsed-data delta-method standard errors from the model-implied multinomial variance of the bin counts, {bf:Var(y_j) ~ y_j}. It reproduces the variance of a regression run in the expanded one-row-per-individual-per-bin data, without building that data, and is the same object as Saez (2010)'s delta method. Efficient and correct when the counts really are multinomial (independent individuals, correct polynomial).

{phang}
{cmd:scale(x2)} multiplies the {cmd:vce(conventional)} covariance by the Pearson overdispersion {bf:phi-hat = sum_j (y_j - yhat_j)^2 / yhat_j / (n_ref - p)} estimated from the reference bins -- the quasi-Poisson standard error, as in {helpb glm}. {cmd:scale(#)} uses a literal factor. Supported only with {cmd:vce(conventional)}; {cmd:phi-hat} is returned in {cmd:e(dispersion)} for every analytic run regardless.

{phang}
{cmd:vce(robust)} (= {cmd:hc1}; {cmd:vce(hc0)}/{cmd:vce(hc2)}/{cmd:vce(hc3)} select the other finite-sample corrections) replaces the multinomial diagonal with the Eicker-White residual meat {bf:(y_j - yhat_j)^2}, HC-corrected. The bunching-mass row, whose fitted residual is ~0 by construction, is instead lifted to {bf:phi-hat * Hstar(1 - Hstar/N)} -- borrowing the reference-region overdispersion into the excluded window (turn this off with {opt nomasscorr}). Robust to arbitrary bin-level heteroskedasticity, to polynomial misspecification treated as noise (round-number heaping, secondary bumps, a neighbouring kink), and to the per-bin part of overdispersion from repeated individuals or year effects. {it:Not} robust to cross-bin correlation (the off-diagonals stay multinomial). Typically wider than {cmd:vce(conventional)} and close to the residual-bootstrap standard errors common in the bunching literature.

{phang}
{cmd:vce(bootstrap)} re-runs the full estimator (integration constraint and all) on {opt bootreps(#)} resampled bin-count vectors and reports the standard deviation of the resulting estimates. The {it:subopts} select the resampling scheme:

{phang2}{cmd:multinomial} (default) -- Dirichlet resample of the bin counts. Equivalent to resampling individuals and re-binning, so it targets the same object as {cmd:vce(conventional)}.{p_end}
{phang2}{cmd:residual} -- resample the reference-bin fit residuals {it:iid} with replacement and add them to the fitted counterfactual (Chetty/CFOP). Imposes a common residual variance across bins, so it is close to {cmd:scale(x2)} in spirit; carries the integration-constraint iteration through each draw.{p_end}
{phang2}{cmd:wild} -- multiply each reference-bin residual by a mean-0 variance-1 weight ({cmd:wildweights(rademacher|mammen|webb)}, default {cmd:rademacher}). Heteroskedasticity-robust: the resampling twin of {cmd:vce(robust)}, and the appropriate version of the residual bootstrap when densities heap.{p_end}
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
Use the bootstrap for inference with 200 repetitions:{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(1) vce(bootstrap) bootreps(200)}{p_end}

{pstd}
Report misspecification-robust (Eicker-White residual) standard errors, comparable to the residual-bootstrap standard errors common in the bunching literature:{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(7) t0(0.2) t1(0.6) vce(robust)}{p_end}

{pstd}
Quasi-Poisson standard errors (one overdispersion parameter):{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(7) t0(0.2) t1(0.6) scale(x2)}{p_end}

{pstd}
Reproduce the Chetty/CFOP residual bootstrap, or its heteroskedasticity-robust (wild) counterpart with percentile intervals:{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(7) t0(0.2) t1(0.6) vce(bootstrap, residual)}{p_end}
{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(7) t0(0.2) t1(0.6) vce(bootstrap, wild percentile) bootreps(999)}{p_end}

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
{synopt:{cmd:e(constant)}}1 if the constant-density response inversion was used, 0 if the exact inversion was used{p_end}
{synopt:{cmd:e(nosplit)}}1 if the bunching mass used the {opt poolmass} calculation, 0 if it used {opt splitmass}{p_end}
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
{synopt:{cmd:e(dispersion)}}Pearson overdispersion {it:phi}-hat from the reference bins, for analytic {cmd:vce()}{p_end}
{synopt:{cmd:e(scalefactor)}}factor actually applied by {opt scale()}, when not 1{p_end}
{synopt:{cmd:e(masscorr)}}1 if the bunching-mass overdispersion correction was applied ({cmd:vce(robust)}/{cmd:hc2}/{cmd:hc3}){p_end}
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
{synopt:{cmd:e(bootci)}}bootstrap CI type: {cmd:normal}, {cmd:bc}, or {cmd:percentile}{p_end}
{synopt:{cmd:e(wildweights)}}wild-bootstrap weight distribution, for {cmd:vce(bootstrap, wild)}{p_end}
{synopt:{cmd:e(scale)}}{opt scale()} specification, when not 1{p_end}
{synopt:{cmd:e(properties)}}usually {cmd:b V} when a variance matrix is posted{p_end}

{synoptset 26 tabbed}{...}
{p2col 5 26 30 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector. With default transformation this contains density and bunching parameters; with {opt notransform} it contains raw estimating-equation coefficients{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of {cmd:e(b)}, when computed{p_end}
{synopt:{cmd:e(table)}}table of binned frequencies and variables used for plotting and diagnostics{p_end}
{synopt:{cmd:e(G)}}delta-method Jacobian for transformed parameters, when available{p_end}
{synopt:{cmd:e(ci_bc)} / {cmd:e(ci_percentile)}}bias-corrected / percentile bootstrap confidence bounds for {cmd:e(b)}, when {cmd:vce(bootstrap, bc)} or {cmd:vce(bootstrap, percentile)} is used{p_end}

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
