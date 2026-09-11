{smcl}
{* *! version date 20260911}{...}
{vieweralsosee "polbunch" "help polbunch"}{...}
{vieweralsosee "polbunchbias" "help polbunchbias"}{...}
{viewerjumpto "Syntax" "polbunch_contrast##syntax"}{...}
{viewerjumpto "Description" "polbunch_contrast##description"}{...}
{viewerjumpto "Options" "polbunch_contrast##options"}{...}
{viewerjumpto "Remarks" "polbunch_contrast##remarks"}{...}
{viewerjumpto "Examples" "polbunch_contrast##examples"}{...}
{viewerjumpto "Stored results" "polbunch_contrast##results"}{...}

{title:Title}

{phang}
{bf:polbunch_contrast} {hline 2} paired binned-bootstrap test of the difference in an estimand between two stored {helpb polbunch} estimates


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:polbunch_contrast} {it:nameA} {it:nameB} [{cmd:,} {it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt target(coef)}}element of {cmd:e(b)} to contrast; default {cmd:bunching:elasticity}{p_end}
{synopt:{opt reps(#)}}bootstrap replications; default {cmd:reps(999)}{p_end}
{synopt:{opt seed(#)}}random-number seed{p_end}
{synopt:{opt boottype(scheme)}}resampling scheme: {cmd:multinomial} (default; {cmd:dirichlet} is a synonym){p_end}
{synopt:{opt level(#)}}confidence level for the reported interval; default {cmd:level(95)}{p_end}
{synopt:{opt nodots}}suppress the replication dots{p_end}
{synoptline}

{p 4 6 2}
{it:nameA} and {it:nameB} are estimates saved with {helpb estimates store} after {cmd:polbunch}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:polbunch_contrast} tests whether a parameter -- by default the implied elasticity -- differs between two {cmd:polbunch} specifications fitted on the {bf:same underlying data}: either the same histogram, or one fitted on a narrower z-range subset of the other's histogram (same bin grid).  It is a postestimation command: the histogram is read from {cmd:e(bins)} of the stored estimates, so nothing needs to be in memory and the data currently loaded (and the active estimates) are left untouched.

{pstd}
The difference {it:d = theta_A - theta_B} is not a simple contrast of two independent numbers: both estimators are computed from the same bin counts and move together.  {cmd:polbunch_contrast} resamples the shared histogram and {bf:re-fits both specifications on every resample}, so the bootstrap distribution of {it:d} carries that dependence directly.  The standard error is the bootstrap SD of {it:d}; because the two estimators are typically strongly correlated (they share the reference-region polynomial fit), that SD is far smaller than the quadrature of the two estimators' own standard errors, and adding the latter would be badly conservative.

{pstd}
This is the general-purpose companion to the {cmd:contrast} option inside {helpb polbunch} (see {help polbunch##contrast:that entry}), which is a shortcut for the common case "my estimator vs. the model-consistent efficient one."  {cmd:polbunch_contrast} compares any two stored {cmd:polbunch} estimates -- different estimators, polynomial orders, inversion or mass axes, or excluded windows.


{marker options}{...}
{title:Options}

{phang}
{opt target(coef)} names the coefficient to contrast.  It must be a column of {cmd:e(b)} in {it:both} stored estimates, e.g. {cmd:bunching:elasticity} (the default), {cmd:bunching:shift}, {cmd:bunching:marginal_response}, {cmd:bunching:number_bunchers}, {cmd:bunching:excess_mass}, {cmd:bunching:delta}.  {cmd:polbunch_contrast} verifies the column exists in both models and errors otherwise (for example {cmd:bunching:delta} is reported by estimator 2 and by estimator 3 under {cmd:constant}/{cmd:poolmass}, but not by estimator 3 on the {cmd:exact}+{cmd:splitmass} axes).

{phang}
{opt reps(#)} sets the number of bootstrap replications (default 999).  Each replication re-fits both specifications, so for estimators 2 and 3 -- which run a delta profile per fit -- this is the dominant cost.

{phang}
{opt seed(#)} seeds the random-number generator.  With the same seed the resamples match those of {cmd:polbunch}'s own {cmd:vce(bootstrap)} / {cmd:contrast}.

{phang}
{opt boottype(scheme)} selects the resampling scheme.  Only {cmd:multinomial} (equivalently {cmd:dirichlet}) is supported -- a Dirichlet reweight of the bin counts, which is the standard binned bootstrap in the bunching literature and equals resampling individuals.  {cmd:residual} and {cmd:wild} resample around a per-model fitted mean and are available only through the {cmd:contrast} option inside {cmd:polbunch}.

{phang}
{opt level(#)} sets the confidence level of the reported percentile interval for the difference.

{phang}
{opt nodots} suppresses the progress dots.


{marker remarks}{...}
{title:Remarks}

{pstd}
{bf:Same data.}  The two estimates must have been fitted on the same histogram, or one on a bin-for-bin subset of the other -- e.g. the same underlying data, but one specification restricted to a narrower z-range (same bin width and alignment).  {cmd:polbunch_contrast} checks {cmd:e(bins)} for this and refuses to proceed otherwise; when the ranges differ it resamples from the wider histogram and refits the narrower specification on only its own bins, so the pairing (and hence the bootstrap SE) still reflects the shared data.  It also notes the range difference in {cmd:r(warning)}, since the difference between the two estimates then partly reflects the choice of window as well as the estimator.  Estimates saved by a {cmd:polbunch} version that predates {cmd:e(bins)} must be re-fitted and re-stored.

{pstd}
{bf:Estimators.}  Any two of estimators 0-4 may be contrasted.  In {cmd:polbunch}'s implementation estimator 4 (Saez) fits a constant per side on the {it:same} reference bins as estimators 1-3, and estimator 0 uses only the left reference bins; the resample feeds each estimator its own bins, so the paired bootstrap covariance is correct in every case.  When the pair mixes counterfactual structures -- the free two-sided fit (0), the restricted model (3), or the Saez two-point trapezoid (4) -- the two elasticities need not share a probability limit even under the model; {cmd:polbunch_contrast} prints a note and sets {cmd:r(mixnote)}.

{pstd}
{bf:Different specifications.}  When the two specifications differ in polynomial order, excluded window or cutoff, {cmd:polbunch_contrast} still runs but prints a note and sets {cmd:r(warning)}: the difference then reflects those choices as well as the estimator.

{pstd}
{bf:Interpretation.}  A rejection says the two specifications give significantly different answers {it:on this sample}; it is not a signed statement that either is biased, and the number is not comparable across other estimator pairs.  Because the two estimators are strongly correlated the test has high power and will flag even economically small gaps -- read the reported difference for magnitude alongside the p-value.  A non-rejection is genuine reassurance that the choice between the two is immaterial in that application.

{pstd}
{bf:Reported tests.}  The normal-approximation p-value uses {it:z = d / SD(d*)}.  The bootstrap p-value is the share of replications with {bf:|d* - d| >= |d|}.  The interval is a percentile interval for {it:d*}.


{marker examples}{...}
{title:Examples}

{pstd}Fit two specifications and store them:{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(7) estimator(2) t0(0.2) t1(0.6) vce(none)}{p_end}
{phang2}{cmd:. estimates store chetty}{p_end}
{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(7) estimator(3) exact splitmass t0(0.2) t1(0.6) vce(none)}{p_end}
{phang2}{cmd:. estimates store efficient}{p_end}

{pstd}Test the elasticity difference, then the shift difference:{p_end}

{phang2}{cmd:. polbunch_contrast chetty efficient, reps(999) seed(20260910)}{p_end}
{phang2}{cmd:. polbunch_contrast chetty efficient, target(bunching:shift)}{p_end}

{pstd}Does the polynomial order matter for the elasticity?{p_end}

{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) polynomial(4) estimator(3) t0(0.2) t1(0.6) vce(none)}{p_end}
{phang2}{cmd:. estimates store deg4}{p_end}
{phang2}{cmd:. polbunch_contrast efficient deg4}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:polbunch_contrast} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(diff)}}{it:nameA} value minus {it:nameB} value of {cmd:target()}{p_end}
{synopt:{cmd:r(value_A)}, {cmd:r(value_B)}}the two stored point values of {cmd:target()}{p_end}
{synopt:{cmd:r(se)}}paired bootstrap standard error of {cmd:r(diff)}{p_end}
{synopt:{cmd:r(corr)}}bootstrap correlation of the two target draws{p_end}
{synopt:{cmd:r(z)}, {cmd:r(p)}}{it:z} = {cmd:r(diff)/r(se)} and its two-sided normal p-value{p_end}
{synopt:{cmd:r(p_bootstrap)}}bootstrap p-value: share of replications with {bf:|d* - d| >= |d|}{p_end}
{synopt:{cmd:r(ci_ll)}, {cmd:r(ci_ul)}}percentile interval for the difference at {cmd:level()}{p_end}
{synopt:{cmd:r(reps)}, {cmd:r(reps_used)}}replications requested and those that yielded a paired difference{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(target)}}the contrasted coefficient{p_end}
{synopt:{cmd:r(name_A)}, {cmd:r(name_B)}}the two stored-estimate names{p_end}
{synopt:{cmd:r(spec_A)}, {cmd:r(spec_B)}}compact labels for the two specifications{p_end}
{synopt:{cmd:r(warning)}}set when the specifications differ in degree / window / cutoff{p_end}
{synopt:{cmd:r(mixnote)}}set when the pair mixes counterfactual structures (estimator 0, 3 or 4){p_end}


{title:Author}

{pstd}Martin Eckhoff Andresen, Department of Economics, University of Oslo.{p_end}
{pstd}martin.eckhoff.andresen@gmail.com{p_end}


{title:Also see}

{psee}
Help:  {helpb polbunch}, {helpb polbunchbias}, {helpb polbunchplot}
{p_end}
