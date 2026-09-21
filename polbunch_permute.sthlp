{smcl}
{* *! version date 20260918}{...}
{vieweralsosee "polbunch" "help polbunch"}{...}
{vieweralsosee "polbunch_contrast" "help polbunch_contrast"}{...}
{vieweralsosee "polbunchbias" "help polbunchbias"}{...}
{viewerjumpto "Syntax" "polbunch_permute##syntax"}{...}
{viewerjumpto "Description" "polbunch_permute##description"}{...}
{viewerjumpto "Options" "polbunch_permute##options"}{...}
{viewerjumpto "Remarks" "polbunch_permute##remarks"}{...}
{viewerjumpto "Examples" "polbunch_permute##examples"}{...}
{viewerjumpto "Stored results" "polbunch_permute##results"}{...}

{title:Title}

{phang}
{bf:polbunch_permute} {hline 2} placebo-cutoff permutation inference for a stored {helpb polbunch} estimate


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:polbunch_permute} [{it:name}] [{cmd:,} {it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt target(coef)}}element of {cmd:e(b)} to permute; default {cmd:bunching:excess_mass}{p_end}
{synopt:{opt step(#)}}cutoff-grid increment, original z units; default the model's own {cmd:e(bw)}{p_end}
{synopt:{opt maxcutoffs(#)}}cap on the total number of candidate placebo cutoffs; default 100, {cmd:0} = unlimited{p_end}
{synopt:{opt deltamax(#)}}upper bound passed to every placebo/observed refit's delta profile; default {cmd:deltamax(1)}{p_end}
{synopt:{opt onesided}}one-sided permutation p-value (default is two-sided on {cmd:|excess_mass|}){p_end}
{synopt:{opt clean}}cross-sample mode: skip the same-data verification, for testing the specification against a different, presumed kink-free sample{p_end}
{synopt:{opt spectests}}also permute whatever Wald/minimum-distance/Hausman/deviance tests the target model already has{p_end}
{synopt:{opt avoid(numlist)}}extra kink locations to keep out of the candidate grid entirely{p_end}
{synopt:{opt level(#)}}percentile spread of the placebo distribution to display; default {cmd:level(95)}{p_end}
{synopt:{opt nodots}}suppress the placebo-refit progress dots{p_end}
{synoptline}

{p 4 6 2}
{it:name} is an estimate saved with {helpb estimates store} after {cmd:polbunch}.  If omitted, {cmd:polbunch_permute} operates on whatever {cmd:polbunch} results are currently active.


{marker description}{...}
{title:Description}

{pstd}
{cmd:polbunch_permute} asks whether the excess mass {cmd:polbunch} found at the true cutoff is unusual, or whether a similarly-sized excess mass would show up at an arbitrary point in the same distribution.  This is a placebo/randomization check in the spirit of Cattaneo, Frandsen and Titiunik (2017, {it:JASA}) local randomization inference for regression discontinuity, and of Ganong-Jager style permutation tests at placebo kink locations: re-fit the identical specification at every viable placebo cutoff across the data, and compare the true cutoff's estimate against that placebo distribution.

{pstd}
The target is fixed to {cmd:bunching:excess_mass} by default and cannot be {cmd:bunching:elasticity}: a placebo location has no real tax or price change, so the elasticity is not defined there, while excess mass -- the gap between the observed and counterfactual mass near the (placebo) cutoff -- is estimable anywhere.

{pstd}
{cmd:polbunch_permute} is a postestimation command, but unlike {helpb polbunch_contrast} it needs the {bf:currently loaded dataset}, not just {cmd:e(bins)}: each placebo cutoff requires its own bin grid (bins are defined relative to the cutoff), so the original {cmd:e(bins)} -- aligned to the true cutoff -- cannot be reused for a different one.  Before doing anything else, it therefore verifies that the data in memory reproduces the stored model's fitted histogram bin-for-bin at the model's own cutoff and bin width, then re-fits directly on the live data for every candidate cutoff.

{pstd}
Every draw -- the true cutoff's and every placebo's -- is {bf:studentized} by its own analytic (delta-method) standard error before comparison.  {cmd:excess_mass}'s sampling variance differs sharply across the income distribution: even after its own density normalization, a placebo cutoff sitting in a sparse region of the data gives a noisier draw than one in a dense region, and comparing raw magnitudes would not account for that.  The permutation p-value and the placebo distribution reported in {cmd:r(placebo_t_*)} are therefore on {cmd:target()/SE(target())}, not raw {cmd:target()} (which is still reported, in {cmd:r(observed)} and {cmd:r(placebo_mean)} etc., for interpretability on its natural scale).


{marker options}{...}
{title:Options}

{phang}
{opt target(coef)} names the coefficient to permute.  Must be a column of {cmd:e(b)} that does not require {cmd:t0()}/{cmd:t1()} -- e.g. {cmd:bunching:excess_mass} (the default), {cmd:bunching:number_bunchers}, {cmd:bunching:shift}, {cmd:bunching:marginal_response}, {cmd:bunching:delta}.  {cmd:target(bunching:elasticity)} is rejected outright (see {help polbunch_permute##description:Description}).

{phang}
{opt step(#)} sets the cutoff-grid increment in original z units.  The default, {cmd:step(0)}, resolves to the stored model's own {cmd:e(bw)}, matching the Ganong-Jager convention of testing at every bin-width-spaced location.

{phang}
{opt maxcutoffs(#)} caps the total number of candidate cutoffs attempted (both directions from the true cutoff combined).  The default is 100; {cmd:maxcutoffs(0)} removes the cap and walks the full usable data range, which can be slow on finely-binned data.

{phang}
{opt deltamax(#)} is passed through to every placebo and observed refit's estimator 2/3 delta-profile search.  It is not reconstructed from the stored model.

{phang}
{opt onesided} switches the permutation p-value from the default two-sided test on {cmd:|target|} to a one-sided test on {cmd:target} itself, matching the theoretical prior that true bunching produces positive excess mass.

{phang}
{opt clean} switches to cross-sample mode: the same-data verification (see {help polbunch_permute##remarks:Remarks}) is skipped, and the specification stored in {it:name} is applied to whatever dataset is currently loaded -- typically a different, presumed kink-free sample, e.g. a pre-policy year or a comparison group with the same variable names.  {cmd:r(observed)} then answers "does spurious bunching show up at this same threshold value in a sample where it shouldn't," a standard placebo check in the kink/bunching literature.  The excluded-window logic still works unchanged (it depends only on {cmd:cutoff}/{cmd:bw}/{cmd:limits}, not on the data), so the nominal cutoff's own donut is still excluded from its placebo comparison set.

{phang}
{opt spectests} permutes whatever specification tests are already sitting in the target model's {cmd:e()} -- nothing new is requested. For each of the omnibus nested-restriction test ({cmd:e(chi2_omnibus)} -- wald for estimator 1/4, minimumdistance for 2/3), Hausman ({cmd:e(chi2_hausman)}), reference-region deviance ({cmd:e(deviance)}), and below-cutoff-only deviance ({cmd:e(deviance_below)}) that is nonmissing on the target, every placebo (and the true-cutoff refit) recomputes that same statistic, and a one-sided permutation p-value ({cmd:r(p_*_perm)}) is reported alongside the target's own asymptotic one. A target fit with an off-label {cmd:test(wald)} on estimator 2/3 (a secondary shape diagnostic, not that estimator's own omnibus test) is permuted separately under {cmd:e(chi2_wald)}/{cmd:r(observed_wald)}/{cmd:r(p_wald_perm)}. See {help polbunch_permute##remarks:Remarks} for the cost and {cmd:norankcheck} caveats.

{phang}
{opt avoid(numlist)} lists extra kink locations -- typically OTHER real kinks elsewhere in the same z-distribution -- that must never enter the candidate grid or any placebo's reference window.  Each point is given the {it:same} excluded-region half-widths as the target model's own donut ({cmd:(cutoff - e(lower_limit))} below, {cmd:(e(upper_limit) - cutoff)} above), recentred at that point -- i.e. "a kink shaped like mine, just centred elsewhere," not a caller-chosen radius.  Filtered at grid construction exactly like the true cutoff's own excluded region: never fit, never counted as attempted/used/failed, tallied instead in {cmd:r(cutoffs_avoid)}.  If a kink's excluded region has a different shape than the target's own, restrict the data directly with this command's own {cmd:if}/{cmd:in} instead (see {help polbunch_permute##remarks:Remarks}).

{phang}
{opt level(#)} sets the percentile spread of the placebo distribution shown in the display and returned in {cmd:r(placebo_lb)}/{cmd:r(placebo_ub)}.

{phang}
{opt nodots} suppresses the progress dots shown while refitting at each placebo cutoff.


{marker remarks}{...}
{title:Remarks}

{pstd}
{bf:Verification.}  By default, {cmd:polbunch_permute} re-bins the current z-variable at the stored model's own cutoff and bin width and requires the result to match {cmd:e(bins)} exactly, row for row.  This confirms the loaded data is the same underlying variable the model was fitted on (possibly finer-grained, e.g. individual-level data behind a pre-binned fit) before any placebo refitting happens.  If the loaded dataset differs -- wrong file, altered observations, a different {cmd:if}/{cmd:in} -- this fails loudly rather than silently returning a placebo distribution built on the wrong data.

{pstd}
{bf:clean (cross-sample mode).}  {opt clean} deliberately disables that check, since its whole point is to apply the stored specification to data the model was {it:not} fitted on.  This trades away the automatic guard against an accidentally-wrong dataset for the flexibility to test a different sample on purpose -- there is no way for {cmd:polbunch_permute} to distinguish "wrong file loaded by mistake" from "correct, deliberately different comparison sample," so the display and {cmd:r(clean)} always flag when this mode was used.  Everything downstream of verification -- the excluded-window/overlap logic, studentization, the grid, the p-value -- is unchanged, because the excluded-region edges are a pure function of {cmd:(cutoff, bw, limits)}, never of the data itself.

{pstd}
{bf:Placebo grid and overlap exclusion.}  Candidate cutoffs step outward from the true cutoff in both directions by {cmd:step()}.  Each candidate is re-fit with the {it:same} number of excluded bins on each side of the cutoff ({cmd:limits(L H)}, recovered from the stored model) and the {it:same} {cmd:allownegative} delta search as the true-cutoff refit.  A candidate is excluded from the reported distribution -- but still counted in {cmd:r(cutoffs_overlap)} -- whenever its own excluded window would overlap the true model's excluded region {cmd:[e(lower_limit), e(upper_limit)]}, so placebo estimates are never contaminated by observations near the real kink.  A candidate that fails to converge, or for which {cmd:target()} is missing, is counted in {cmd:r(cutoffs_skipped)} instead.

{pstd}
{bf:Other real kinks: avoid() vs. if/in.}  If a neighbouring kink's excluded region is shaped like the target's own (same donut half-widths), just list its location in {opt avoid()} -- candidates near it are filtered out cheaply, before ever being fit, and counted in {cmd:r(cutoffs_avoid)} ({cmd:r(avoid_user)}).  If it isn't the same shape, exclude its data directly with this command's own {cmd:if}/{cmd:in} instead.  Either way, a candidate whose window merely {it:touches} the resulting hole is excluded outright, {it:not} just when the hole happens to break that specific candidate's fit: before the candidate grid is built, {cmd:polbunch_permute} scans the sorted loaded data (one cheap pass, no model fits) for any maximal run of CONSECUTIVE observations that all fail your {cmd:if}/{cmd:in} -- exact edges, independent of this model's own bin width or phase (deliberately not a bin-count scan: two kinks need not sit on the same bw-spaced phase, and a phased bin scan would under-detect the hole by up to a full bin width right at the boundary when they don't) -- and folds each such run into the identical pre-fit window-overlap filter {cmd:avoid()} uses -- counted in {cmd:r(cutoffs_avoid)} ({cmd:r(avoid_auto)}), never fit, never counted as attempted/used/failed.  This matters because a candidate whose window straddles the hole but still converges on whatever data survives outside it is {it:not} a clean placebo location -- its own reference-region fit is missing a real chunk of data, exactly the contamination the excluded-region logic exists to rule out (see the overlap exclusion above), even though it didn't happen to crash.  Data that was already absent regardless of your {cmd:if}/{cmd:in} is left alone (ordinary sparse-data tolerance is unaffected).  Anything that still fails to fit and isn't explained by an auto-detected gap falls back to a cheaper diagnostic: the candidate is re-fit once more on the identical window with your {cmd:if}/{cmd:in} lifted, and if that refit converges, the original failure is attributed to your restriction rather than to the placebo location -- dropped from {cmd:r(placebo_draws)} entirely and tallied in {cmd:r(cutoffs_candif)} instead (its fitted value is discarded either way, since it would still be contaminated by whatever your {cmd:if}/{cmd:in} was excluding).  A candidate that fails even without your restriction is a genuine non-convergence and is still counted in {cmd:r(cutoffs_skipped)}.

{pstd}
{bf:allownegative.}  The true-cutoff specification is re-fit once too (to produce {cmd:r(observed)}), with {cmd:allownegative} forced on just like every placebo -- a placebo location has no structural reason for the delta &gt;= 0 constraint that {cmd:polbunch} otherwise imposes by default, so imposing it only at the true cutoff would make {cmd:r(observed)} incomparable to the placebo draws.

{pstd}
{bf:Determinism.}  This is a deterministic grid walk, not a resampling procedure -- there is no {cmd:seed()} and no randomness.  The same call on the same data always returns the same {cmd:r(placebo_draws)}.

{pstd}
{bf:Permutation p-value.}  {cmd:r(p)} uses the standard "+1" exact-permutation correction: the count of studentized placebo draws at least as extreme as {cmd:r(observed_t)}, plus one, divided by the number of used placebo draws, plus one -- this avoids reporting an implausible {cmd:p=0} off a finite candidate grid.

{pstd}
{bf:Studentization.}  Every refit uses {cmd:vce(conventional)} (a cheap analytic/delta-method SE), independent of whatever {cmd:vce()} the original stored model used -- a bootstrap or cluster {cmd:vce()} would be far too costly to repeat once per placebo cutoff.  A candidate whose SE for {cmd:target()} is missing or non-positive is treated as a failed refit ({cmd:r(cutoffs_skipped)}), the same as a candidate that fails to converge.

{pstd}
{bf:Scope limits.}  {cmd:nodrop}, {cmd:nozero}, {cmd:positive}, {cmd:scale()}, {cmd:nomasscorr} and {cmd:vce(cluster ...)} are not reconstructed from the stored model (none of these survive to the final {cmd:e()} of a {cmd:polbunch} fit).  If the original fit used {cmd:nodrop} or {cmd:nozero} in a way that changes the bin grid itself, the verification step above will fail rather than silently produce a mismatched placebo distribution.

{pstd}
{bf:spectests.}  Detection reads {cmd:e()} directly ("what's there"), not the target's typed {cmd:test()} option -- if a Wald/minimum-distance/Hausman chi2 is present, every placebo (and the observed refit) is re-fit with whichever {cmd:test()} value covers exactly the union of what was found ({cmd:test(all)} when more than one type is present); deviance and below-cutoff deviance need no {cmd:test()} at all, since {cmd:vce(conventional)} already posts them on every fit. The p-values are always one-sided (chi2 and deviance are non-negative; larger is more extreme), regardless of {cmd:onesided}. The target's own null-imposed bootstrap CDF ({cmd:test(..., boot)}) is never re-run per placebo -- only the target model computes that, once.

{pstd}
{bf:spectests cost and the norankcheck conflict.}  Deviance/below-cutoff deviance are free -- already posted by every {cmd:vce(conventional)} fit. Wald/minimum-distance/Hausman are not: {cmd:test()} triggers an internal restricted-vs-unrestricted refit inside every placebo's {cmd:polbunch} call, roughly doubling its cost, so {cmd:maxcutoffs()} is worth tightening when any chi2 test is requested. More importantly, {cmd:polbunch} itself forces {cmd:test(none)} internally whenever {cmd:norankcheck} is given (documented in {cmd:polbunch.ado}: the tests need an identified unrestricted comparison fit, which {cmd:norankcheck} deliberately forgoes) -- so whenever {cmd:spectests} detects a chi2 test, {cmd:norankcheck} is dropped from every placebo/observed refit for that call. This reopens exactly the polynomial-degree-drift risk {cmd:norankcheck} exists to close (see {help polbunch_contrast:polbunch_contrast}'s use of the same option): a placebo cutoff could, in principle, come back at a lower degree than the target's own. {cmd:r(spectest_degree_drift)} counts how many used placebo cutoffs actually did so -- 0 is the expected/common case, but it is reported rather than assumed. Deviance-only requests (no chi2 test present) are unaffected and keep {cmd:norankcheck}.


{marker examples}{...}
{title:Examples}

{pstd}Fit and store a model, then test its excess mass against placebo cutoffs:{p_end}

{phang2}{cmd:. polbunch z, cutoff(50000) bw(500) polynomial(7) estimator(3) exact splitmass vce(none)}{p_end}
{phang2}{cmd:. estimates store main}{p_end}
{phang2}{cmd:. polbunch_permute main}{p_end}

{pstd}One-sided test, coarser grid, on the currently active (unstored) estimate:{p_end}

{phang2}{cmd:. polbunch_permute, onesided step(1000) maxcutoffs(40)}{p_end}

{pstd}Placebo check against a pre-policy (presumed kink-free) year, same variable names:{p_end}

{phang2}{cmd:. use policy_year, clear}{p_end}
{phang2}{cmd:. polbunch z, cutoff(50000) bw(500) polynomial(7) estimator(3) exact splitmass vce(none)}{p_end}
{phang2}{cmd:. estimates store main}{p_end}
{phang2}{cmd:. use preperiod_year, clear}{p_end}
{phang2}{cmd:. polbunch_permute main, clean}{p_end}

{pstd}Also permute the target's own specification tests:{p_end}

{phang2}{cmd:. polbunch z, cutoff(50000) bw(500) polynomial(7) estimator(3) exact splitmass vce(conventional) test(omnibus)}{p_end}
{phang2}{cmd:. estimates store main}{p_end}
{phang2}{cmd:. polbunch_permute main, spectests}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:polbunch_permute} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(observed)}}{cmd:target()} at the true cutoff, refit with {cmd:allownegative}{p_end}
{synopt:{cmd:r(observed_se)}}analytic SE of {cmd:r(observed)}{p_end}
{synopt:{cmd:r(observed_t)}}{cmd:r(observed)} / {cmd:r(observed_se)} -- the statistic actually permuted{p_end}
{synopt:{cmd:r(p)}}permutation p-value on the studentized statistic (one- or two-sided per {cmd:onesided}){p_end}
{synopt:{cmd:r(placebo_mean)}, {cmd:r(placebo_sd)}}mean / sd of the used placebo draws, raw scale{p_end}
{synopt:{cmd:r(placebo_min)}, {cmd:r(placebo_max)}}range of the used placebo draws, raw scale{p_end}
{synopt:{cmd:r(placebo_lb)}, {cmd:r(placebo_ub)}}{cmd:level()}% percentile range of the used placebo draws, raw scale{p_end}
{synopt:{cmd:r(placebo_t_mean)}, {cmd:r(placebo_t_sd)}}mean / sd of the used placebo draws, studentized{p_end}
{synopt:{cmd:r(placebo_t_lb)}, {cmd:r(placebo_t_ub)}}{cmd:level()}% percentile range of the used placebo draws, studentized{p_end}
{synopt:{cmd:r(n_placebo)}}number of placebo cutoffs used (converged, usable SE, no overlap){p_end}
{synopt:{cmd:r(n_placebo_tried)}}total number of candidate placebo cutoffs attempted{p_end}
{synopt:{cmd:r(cutoffs_skipped)}}candidates that genuinely failed to converge or lacked {cmd:target()}/a usable SE{p_end}
{synopt:{cmd:r(cutoffs_overlap)}}candidates excluded for overlapping the true excluded region{p_end}
{synopt:{cmd:r(cutoffs_avoid)}}candidates excluded for overlapping an {cmd:avoid()} region or an auto-detected {cmd:if}/{cmd:in} gap{p_end}
{synopt:{cmd:r(avoid_user)}}number of explicit {cmd:avoid()} points{p_end}
{synopt:{cmd:r(avoid_auto)}}number of contiguous data gaps auto-detected from the caller's {cmd:if}/{cmd:in}{p_end}
{synopt:{cmd:r(cutoffs_candif)}}candidates excluded because their window only failed to fit inside the caller's own {cmd:if}/{cmd:in} (not explained by an auto-detected gap){p_end}
{synopt:{cmd:r(step)}}the cutoff-grid increment actually used{p_end}
{synopt:{cmd:r(bw)}}the stored model's bin width{p_end}
{synopt:{cmd:r(cutoff_true)}}the true cutoff{p_end}
{synopt:{cmd:r(onesided)}}1 if {cmd:onesided} was specified{p_end}
{synopt:{cmd:r(clean)}}1 if {cmd:clean} was specified (same-data verification skipped){p_end}
{synopt:{cmd:r(level)}}the {cmd:level()} used{p_end}
{synopt:{cmd:r(spectests)}}1 if {cmd:spectests} was specified{p_end}
{synopt:{cmd:r(observed_omnibus)}, {cmd:r(p_omnibus_perm)}}the target's own omnibus (nested-restriction) chi2 at the true cutoff and its permutation p-value; only when {cmd:e(chi2_omnibus)} was present on the target{p_end}
{synopt:{cmd:r(observed_wald)}, {cmd:r(p_wald_perm)}}Wald chi2 and its permutation p-value; only when the target used an off-label {cmd:test(wald)} on estimator 2/3, so {cmd:e(chi2_wald)} was present (never both this and {cmd:r(observed_omnibus)} for the same target){p_end}
{synopt:{cmd:r(observed_hausman)}, {cmd:r(p_hausman_perm)}}Hausman chi2 and its permutation p-value; only when {cmd:e(chi2_hausman)} was present{p_end}
{synopt:{cmd:r(observed_deviance)}, {cmd:r(p_deviance_perm)}}reference-region deviance and its permutation p-value; only when {cmd:e(deviance)} was present{p_end}
{synopt:{cmd:r(observed_deviance_below)}, {cmd:r(p_deviance_below_perm)}}below-cutoff-only deviance and its permutation p-value; only when {cmd:e(deviance_below)} was present{p_end}
{synopt:{cmd:r(spectest_degree_drift)}}with {cmd:spectests}: count of used placebo cutoffs whose polynomial degree differs from the target's (see Remarks){p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(target)}}the permuted coefficient{p_end}
{synopt:{cmd:r(name)}}the stored-estimate name used (empty if the active estimate was used){p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(placebo_draws)}}one row per candidate cutoff: cutoff, {cmd:target()} value (missing if unused), SE (missing if unused), status (1 = used, 0 = failed, -1 = overlap-filtered){p_end}
{synopt:{cmd:r(spectest_draws)}}with {cmd:spectests}: one row per candidate cutoff -- cutoff, chi2_omnibus, chi2_wald, chi2_hausman, deviance, deviance_below, polynomial degree used (missing where not applicable or not used){p_end}


{title:Author}

{pstd}Martin Eckhoff Andresen, Department of Economics, University of Oslo.{p_end}
{pstd}martin.eckhoff.andresen@gmail.com{p_end}


{title:Also see}

{psee}
Help:  {helpb polbunch}, {helpb polbunch_contrast}, {helpb polbunchbias}, {helpb polbunchplot}
{p_end}
