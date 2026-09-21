*! polbunch_permute version date 20260918
* Author: Martin Eckhoff Andresen
* This program is part of the polbunch package.
*
* Postestimation: placebo-cutoff permutation/randomization inference for a
* stored (or currently active) polbunch estimate, in the spirit of
* Cattaneo-Frandsen-Titiunik (2017, JASA) local randomization inference and
* Ganong-Jager style permutation tests at placebo kink locations.
*
* By default verifies the data currently in memory reproduces the stored
* model's histogram bin-for-bin at its own cutoff/bw (same-sample mode),
* then re-fits the identical specification -- with allownegative forced
* on, no tax rates (no real tax/price change exists at a placebo
* location), and a cheap analytic vce() so every draw can be studentized
* by its own delta-method SE (excess_mass's sampling variance differs
* sharply across the income distribution even after its own density
* normalization) -- at every viable placebo cutoff across the data,
* skipping any candidate whose own excluded window would overlap the
* true model's excluded region.  Reports a permutation p-value for the
* true cutoff's studentized bunching:excess_mass (or another target())
* against the studentized placebo distribution.
*
* Window recentering: every PLACEBO cutoff's estimation window is RECENTERED
* on its own cutoff, using the left/right halfwidths e(zlo)/e(zhi) say the
* target model actually reached -- the outer edges of Hstored=e(bins),
* i.e. what the fit truly consumed, however if/in/nodrop/nozero produced it,
* and independent of whatever happens to be loaded now. This replaces
* replaying the target's literal if/in bounds unchanged at every placebo
* (which silently starves off-centre placebos of data on one side, caps the
* whole candidate grid at the target's own literal window even when the
* loaded data reaches further, and can leave a candidate that overlaps the
* true excluded region undetected as such only by accident). The
* TRUE-CUTOFF refit is untouched -- it keeps if/in literally, since it
* doubles as the bin-for-bin verification against Hstored and must
* reproduce the target's exact original sample, not a bin-extent
* approximation of it. Tradeoff, accepted for simplicity: if if/in also
* carried a substantive sample filter beyond windowing (e.g. an occupation
* restriction), that filter is not replayed at placebo cutoffs -- only the
* window width is. Works identically for individual-level and pre-binned
* (freq zvar) data, since e(zlo)/e(zhi) are set the same way for both.
*
* clean: cross-sample mode for testing a specification against a
* deliberately DIFFERENT, presumed-kink-free sample (e.g. a pre-policy
* year) -- skips the bin-for-bin same-data verification (the whole
* point is that the loaded data is NOT the data the model was fitted
* on), while everything else -- including the self-window exclusion
* around the nominal cutoff, which is a pure function of
* (cutoff,bw,limits) and so needs no change -- works unmodified.
*
* spectests: permutes whatever specification tests the TARGET model
* already has posted in e() -- its own omnibus (wald for estimator 1/4,
* minimumdistance for 2/3, e(chi2_omnibus)) / Hausman chi2 (from its own
* test() choice) and the reference-region deviance / below-cutoff-only
* deviance (posted by every polbunch fit that isn't vce(none)/
* multinomial|bayesian-boot, i.e. free -- see polbunch_estat). A target
* fit with an off-label test(wald) on estimator 2/3 (a secondary shape
* diagnostic, NOT that estimator's own omnibus test -- see polbunch's
* test() docs) is detected and permuted separately, under its own
* e(chi2_wald), never folded into the omnibus slot. Nothing new is
* computed; each detected field gets its own placebo distribution and a
* one-sided permutation p-value (chi2/deviance are non-negative, larger
* = more extreme, so "onesided" doesn't apply here). The chi2 tests are
* NOT free: test() triggers an internal restricted-vs-unrestricted refit
* inside every placebo's polbunch call, roughly doubling its cost;
* deviance is free (already posted on every conventional-vce fit). The
* target model's own null-imposed bootstrap CDF (test(..., boot)) is
* never re-run per placebo -- a bootstrap inside a permutation loop
* defeats the point of both.
*
*   polbunch_permute [name] [if] [in] [, target(bunching:excess_mass) step(#) ///
*                               maxcutoffs(100) deltamax(1) onesided ///
*                               clean spectests avoid(numlist) level(#) nodots ]
*
* if/in: restricts which of the CURRENTLY LOADED data may be used when
* building the placebo-cutoff grid and fitting each placebo -- ANDed onto
* the data-range summarize that determines how far candidates may extend
* (zdlo/zdhi) and onto every placebo's own recentered window
* (inrange(zvar,lo_i,hi_i)). It does NOT touch the TRUE-cutoff refit, which
* must keep using exactly inrange(zvar,e(zlo),e(zhi)) to verify the loaded
* data reproduces the target's histogram bin-for-bin -- if this if/in would
* itself change that histogram (e.g. it restricts away part of the sample
* the target was fitted on), pass clean to skip that check. If if/in empties
* one or more CONTIGUOUS stretches of the raw data outright (e.g.
* !inrange(zvar,lo,hi) carving out another kink), those stretches are
* auto-detected (one cheap data scan on sorted z -- "gap spells," exact and
* independent of this model's own bw-phased bin grid, not a bw-width bin
* scan -- see the code for why phase independence matters when two kinks
* don't share the same bin phase) and folded into the SAME pre-fit
* window-overlap filter as avoid() -- a candidate whose window merely
* touches the hole is excluded before ever being fit, not just when the
* hole happens to break that candidate's specific fit. This matters: a
* candidate whose window touches the hole but still converges on its
* surviving data is NOT the same as a clean placebo location -- see avoid()
* below for why. Only genuinely if/in-created gaps trigger this (data
* already absent regardless of if/in is left alone); anything else that
* still fails to fit falls back to the reactive candif classification
* described further down.
*
* avoid(#1 #2 ...): extra kink locations (e.g. OTHER real kinks in the same
* z-distribution) that must not enter the candidate grid or any placebo's
* reference window -- the whole point of if/in-excluding them from the data
* in the first place. Each point is given the SAME excluded-region
* half-widths as the target model's own donut ((cutoff-lower_limit),
* (upper_limit-cutoff)) recentred at that point -- i.e. "a kink shaped like
* mine, just centred elsewhere" -- not a caller-specified radius; pass if/in
* instead for a differently-shaped exclusion (its own gap is then
* auto-detected and filtered exactly the same way, see above). Filtered at
* grid construction, exactly like the true cutoff's own excluded region --
* never fit, never counted as attempted/used/failed, tallied in
* r(cutoffs_avoid) instead (r(avoid_user)/r(avoid_auto) split the count by
* explicit avoid() points vs. auto-detected if/in gaps).
*
* Returns r(p), r(observed), r(observed_se), r(observed_t), r(n_placebo),
* r(n_placebo_tried), r(cutoffs_skipped), r(cutoffs_overlap),
* r(cutoffs_avoid), r(avoid_user), r(avoid_auto), r(cutoffs_candif),
* r(placebo_mean/_sd/_min/_max/_lb/_ub) (raw scale),
* r(placebo_t_mean/_sd/_lb/_ub) (studentized scale), r(step), r(bw),
* r(cutoff_true), r(onesided), r(clean), r(spectests), r(level),
* r(window_recentered), r(window_leftdist), r(window_rightdist),
* r(target), r(name), r(placebo_draws) (matrix, named columns cutoff,
* value, se, status [1=used, 0=failed]; sorted ascending by cutoff).
* Candidates whose own recentered window would overlap the true excluded
* region are excluded before this matrix is even built (counted in
* r(cutoffs_overlap) instead) -- n_placebo_tried and this matrix's row
* count both already exclude them, so there is no -1/overlap status value
* to see here (status -1 remains theoretically possible only in the
* window_recentered==0 fallback, where overlap can only be discovered
* after fitting). Candidates that fail to fit ONLY because of the
* caller's own if/in (e.g. carving out ANOTHER real kink's excluded
* region) are diagnosed by a cheap refit of the same window without that
* if/in; if that refit converges, the candidate is dropped from
* placebo_draws entirely (never counted as attempted, used, or failed)
* and tallied in r(cutoffs_candif) instead -- so a wide if/in-driven data
* gap elsewhere in the sample is never misreported as non-convergence.
* With spectests: for each of omnibus/wald/hausman/deviance/
* deviance_below present on the target, r(observed_<name>) and
* r(p_<name>_perm), plus matrix r(spectest_draws) (named columns cutoff,
* chi2_omnibus, chi2_wald, chi2_hausman, deviance, deviance_below,
* polynomial; same row order as placebo_draws -- missing where not
* applicable). chi2_omnibus is the target's own nested-restriction test
* (wald for estimator 1/4, minimumdistance for 2/3); chi2_wald here is
* ONLY the off-label diagnostic case (test(wald) explicitly requested on
* estimator 2/3), never the same thing as chi2_omnibus for the same fit.

cap program drop polbunch_permute
program define polbunch_permute, rclass
	version 16.0

	syntax [anything] [if] [in] [, ///
		TARGET(string) ///
		STEp(real 0) ///
		MAXCutoffs(integer 100) ///
		DELTAmax(real 1) ///
		ONESided ///
		CLEAN ///
		SPECtests ///
		AVOID(numlist) ///
		Level(cilevel) ///
		NODOTS ]

	/* snapshot this command's own if/in immediately -- the second syntax   *
	 * call below (recovering the ORIGINAL model's if/in from e(cmdline))   *
	 * reuses and overwrites the `if'/`in' locals, so these must be copied  *
	 * out first or they'd be silently lost.                                */
	local candif `"`if'"'
	local candin `"`in'"'
	local candifexpr ""
	if `"`candif'"' != "" local candifexpr = subinstr(`"`candif'"', "if ", "", 1)

	local name = strtrim(`"`anything'"')
	local nname : word count `name'
	if `nname' > 1 {
		di as error "specify at most one stored estimate name."
		exit 198
	}

	/* target -- default excess_mass; elasticity is categorically rejected
	   since t0()/t1() are deliberately omitted from every refit here (no
	   real tax/price change exists at a placebo location). */
	if "`target'" == "" local target "bunching:excess_mass"
	local target = strtrim("`target'")
	if "`target'" == "bunching:elasticity" {
		di as error "target(bunching:elasticity) is not supported here: elasticity requires"
		di as error "t0()/t1(), which are deliberately omitted from every placebo/true-cutoff"
		di as error "refit (there is no real tax or price change at a placebo location).  Use"
		di as error "the default target(bunching:excess_mass), or bunching:number_bunchers,"
		di as error "bunching:shift, bunching:marginal_response, bunching:delta."
		exit 198
	}
	local tlab : subinstr local target "bunching:" "", all
	local tlab : subinstr local tlab "_" " ", all

	if `step' < 0 {
		di as error "step() must be nonnegative (0 = use the model's own bin width)."
		exit 198
	}
	if `maxcutoffs' < 0 {
		di as error "maxcutoffs() must be nonnegative (0 = unlimited)."
		exit 198
	}
	if `deltamax' <= 0 {
		di as error "deltamax() must be a positive number."
		exit 198
	}

	/* ------------------------------------------------------------------ *
	 *  Locate the model to permute, and set up how to restore whatever   *
	 *  was active before this command on exit.  When a name is given,    *
	 *  the currently active estimates are unrelated -- _estimates hold   *
	 *  stashes them away, then estimates restore loads the named one.    *
	 *  When no name is given, the target model IS the currently active   *
	 *  one, so it must stay live in e() for the extraction below; it is  *
	 *  snapshotted with estimates store instead (hold would immediately  *
	 *  clear e() out from under us) and restored the same way at the end.*
	 * ------------------------------------------------------------------ */
	tempname _cur
	local curkind "hold"
	if "`name'" != "" {
		_estimates hold `_cur', restore nullok
		capture estimates restore `name'
		if _rc {
			di as error "could not restore `name' -- use {cmd:estimates store `name'} after fitting."
			exit 198
		}
	}
	else {
		if "`e(cmd)'" != "polbunch" {
			di as error "no polbunch estimation results are currently active; fit a model first, or give a stored name."
			exit 301
		}
		quietly estimates store `_cur'
		local curkind "store"
	}
	if "`e(cmd)'" != "polbunch" {
		if "`name'" != "" di as error "`name' is not a polbunch estimate."
		else di as error "no polbunch estimation results are currently active; fit a model first, or give a stored name."
		exit 301
	}
	if !inlist(e(estimator),0,1,2,3,4) {
		di as error "the stored/active estimate has an unrecognised e(estimator)."
		exit 301
	}
	capture confirm matrix e(bins)
	if _rc {
		di as error "this estimate has no e(bins); it was fitted with an older polbunch."
		di as error "Re-fit it (and, if using a stored name, {cmd:estimates store} it again)."
		exit 198
	}

	/* ------------------------------------------------------------------ *
	 *  spectests: detect which test fields the TARGET already has --     *
	 *  nothing new requested, just "what's there." Read before anything  *
	 *  else touches e().                                                  *
	 * ------------------------------------------------------------------ */
	/* hadomni/hadwald are mutually exclusive for the same target fit: a
	   normal fit's own nested-restriction test posts e(chi2_omnibus)
	   (whichever of wald/minimumdistance mechanically applies); only an
	   explicit off-label test(wald) on estimator 2/3 posts e(chi2_wald)
	   instead, and it never runs alongside hausman (see polbunch's
	   test() option) -- so hadwald & hadhaus never co-occur either. */
	local hadomni = 0
	local hadwald = 0
	local hadhaus = 0
	local haddev  = 0
	local haddevb = 0
	if "`spectests'" != "" {
		local hadomni = !missing(e(chi2_omnibus))
		local hadwald = !missing(e(chi2_wald))
		local hadhaus = !missing(e(chi2_hausman))
		local haddev  = !missing(e(deviance))
		local haddevb = !missing(e(deviance_below))
	}
	local nchi2types = `hadomni' + `hadwald' + `hadhaus'
	if `hadomni' & `hadhaus' local spectestopt "test(all)"
	else if `hadomni' local spectestopt "test(omnibus)"
	else if `hadwald' local spectestopt "test(wald)"
	else if `hadhaus' local spectestopt "test(hausman)"
	else local spectestopt "test(none)"

	/* polbunch itself forces test(none) internally whenever norankcheck   *
	 * is given (a documented design choice: Wald/MD/Hausman need an       *
	 * identified unrestricted comparison fit, which norankcheck forgoes   *
	 * by construction) -- so norankcheck and a requested chi2 test are    *
	 * mutually exclusive. Deviance/deviance_below don't go through        *
	 * test() at all and are unaffected, so norankcheck only needs to      *
	 * drop when a chi2 test was actually detected. Dropping it trades     *
	 * away the "same K at every placebo" guarantee for those placebos     *
	 * only -- e(polynomial) is tracked below so degree drift is visible,  *
	 * not silent. */
	local rankopt "norankcheck"
	if `nchi2types' > 0 local rankopt ""

	/* ------------------------------------------------------------------ *
	 *  Extract everything needed from e() before anything else touches   *
	 *  it -- estimator/spec toggles, the fitted histogram, and the       *
	 *  excluded-region edges (with a fallback to estimator 4's own       *
	 *  scalar names, though the main names are posted for every          *
	 *  estimator).                                                       *
	 * ------------------------------------------------------------------ */
	tempname Hstored
	matrix `Hstored' = e(bins)

	local estimator0 = e(estimator)
	local k0         = e(polynomial)
	local cut0       = e(cutoff_orig)
	local bw0        = e(bw)
	local con0       = e(constant)
	local nsp0       = e(nosplit)
	local islog0     = (e(log) == 1)
	local nonorm0    = ("`e(normalize)'" == "nonormalize")

	local lo0 = e(lower_limit)
	local hi0 = e(upper_limit)
	if missing(`lo0') | missing(`hi0') {
		local lo0 = e(zL_excl_orig)
		local hi0 = e(zH_excl_orig)
	}
	if missing(`lo0') | missing(`hi0') {
		di as error "could not recover the excluded-region edges (e(lower_limit)/e(upper_limit))"
		di as error "from this estimate."
		exit 198
	}

	/* avoid(): extra kink locations to keep out of the candidate grid,     *
	 * beyond the target's own (lo0,hi0]. Each point gets the SAME excluded- *
	 * region half-widths as the target's own donut -- i.e. it is treated   *
	 * as "a kink shaped like mine, just centred elsewhere" -- not a         *
	 * caller-specified radius. That is a real assumption: if another kink's *
	 * own excluded region is a different width, use if/in on this command  *
	 * instead (see the note on candif-classified candidates below), which  *
	 * accepts an arbitrary caller-defined exclusion instead of this shape-  *
	 * copying shortcut.                                                     */
	local navoid = 0
	if `"`avoid'"' != "" {
		local donutlo0 = `cut0' - `lo0'
		local donuthi0 = `hi0' - `cut0'
		foreach p of numlist `avoid' {
			local ++navoid
			local avoidlo`navoid' = `p' - `donutlo0'
			local avoidhi`navoid' = `p' + `donuthi0'
		}
	}
	local navoiduser = `navoid'

	/* recover limits(L H) in bin units from the stored window edges,      *
	 * exactly mirroring polbunch_contrast's edge-vs-crossbin detection.   */
	local m0 = `Hstored'[1,2]
	local frac = mod((`cut0' - `m0')/`bw0', 1)
	if `frac' > 0.5 local frac = `frac' - 1
	local edge0 = (abs(abs(`frac') - 0.5) < 0.25)
	if `edge0' {
		local L0    = round((`cut0' - `lo0')/`bw0')
		local Hlim0 = round((`hi0' - `cut0')/`bw0')
	}
	else {
		local L0    = round((`cut0' - `lo0')/`bw0' - 0.5)
		local Hlim0 = round((`hi0' - `cut0')/`bw0' - 0.5)
	}

	/* recover the original varlist / if / in by re-parsing e(cmdline) --  *
	 * the outer polbunch program's own final ereturn post always resets   *
	 * e(depvar) to "freq", so this is the only way to get them back.      */
	local ecmdline `"`e(cmdline)'"'
	gettoken _junk _rest : ecmdline
	local 0 `"`=strtrim(`"`_rest'"')'"'
	capture syntax varlist(min=1 max=2) [if] [in], [*]
	if _rc {
		di as error "could not recover the original varlist/if/in from e(cmdline):"
		di as error `"  `ecmdline'"'
		di as error "This estimate may predate this feature; re-fit it with the current polbunch."
		exit 198
	}
	local zvarlist `"`varlist'"'
	local zif `"`if'"'
	local zin `"`in'"'
	local zifexpr ""
	if `"`zif'"' != "" local zifexpr = subinstr(`"`zif'"', "if ", "", 1)
	foreach v of local zvarlist {
		capture confirm variable `v'
		if _rc {
			di as error "variable `v' (from the original polbunch specification) is not found"
			di as error "in the currently loaded dataset."
			exit 111
		}
	}
	local nvars0 : word count `zvarlist'
	local bwopt = cond(`nvars0'==1, "bw(`bw0')", "")
	local zvar : word `nvars0' of `zvarlist'

	/* ------------------------------------------------------------------ *
	 *  Window bounds, ALWAYS read from e(zlo)/e(zhi) -- the actual outer  *
	 *  edges of the bins in Hstored=e(bins) (already bin-midpoint +/- a   *
	 *  half bin-width, i.e. true edges, not midpoints -- verified: refit  *
	 *  with inrange(zvar,e(zlo),e(zhi)) reproduces Hstored to machine     *
	 *  precision), i.e. exactly what the target model fit, regardless of  *
	 *  how if/in/nodrop/nozero produced it and regardless of what happens *
	 *  to be loaded now (extra rows appended, a different dataset under   *
	 *  the same varnames, ...). if/in is NOT replayed anywhere below --   *
	 *  every refit (the true cutoff AND every placebo) uses ONLY          *
	 *  inrange(zvar, lo, hi) for its window, `nodrop' forced in           *
	 *  `basespec' so that boundary doesn't get further, recursively       *
	 *  narrowed by polbunch's own edge-trim (which operates on THIS       *
	 *  call's raw min/max, not the target's -- reapplying if/in's own     *
	 *  bounds as a fresh restriction on continuous data would otherwise   *
	 *  shave off an extra partial bin on each side, every time).          *
	 *  Tradeoff, accepted for simplicity: if if/in carried a substantive  *
	 *  sample filter (e.g. an occupation restriction) beyond windowing,   *
	 *  that filter is NOT replayed anywhere -- only the window e(zlo)/    *
	 *  e(zhi) records is. If the currently loaded data has no more room   *
	 *  than the target's own window (e.g. the target already used its    *
	 *  entire dataset), there is nowhere to place a placebo at all -- the *
	 *  existing "no placebo cutoffs fit" error below fires, by design.    */
	local zlo0 = e(zlo)
	local zhi0 = e(zhi)
	local recenter = !missing(`zlo0') & !missing(`zhi0')
	local leftdist0  = `cut0' - `zlo0'
	local rightdist0 = `zhi0' - `cut0'

	local inv0  = cond(`con0',"constant","exact")
	local mass0 = cond(`nsp0',"poolmass","splitmass")
	local lg0   = cond(`islog0',"log","")
	local nn0   = cond(`nonorm0',"nonormalize","")

	/* nodrop forced on every reconstructed call (true cutoff AND placebos)  *
	 * -- NOT because the target itself necessarily used nodrop, but        *
	 * because every reconstructed window below is expressed as an exact    *
	 * inrange(zvar,lo,hi) edge pair (from e(zlo)/e(zhi) at the true         *
	 * cutoff, translated at each placebo) and polbunch's own edge-trim     *
	 * operates on THIS call's raw min/max, not the target's -- without     *
	 * nodrop it would shave another partial bin off whichever edge isn't   *
	 * exactly hit by a sampled value, silently narrowing an already-exact  *
	 * window (see the note above). nozero/positive/scale()/nomasscorr/     *
	 * vce(cluster ...) are still not reconstructed (none survive to final  *
	 * e()); if the target used those in a way that changes the bin grid    *
	 * itself, the verification refit below fails loudly rather than        *
	 * silently mismatching.                                                *
	 *                                                                      *
	 * vce(conventional): excess_mass's sampling variance differs sharply  *
	 * across placebo locations (sparser regions of the distribution give  *
	 * noisier draws even after excess_mass's own density normalization),  *
	 * so every draw is studentized by its own delta-method SE before      *
	 * comparison -- a cheap analytic SE, not whatever (possibly           *
	 * bootstrap) vce() the headline call used, since that cost would be   *
	 * paid once per placebo cutoff.                                       *
	 *                                                                      *
	 * test(): `spectestopt' is test(none) unless spectests found chi2      *
	 * fields on the target, in which case it reproduces exactly the       *
	 * union of what's there (test(all) covers omnibus+hausman; the        *
	 * single omnibus/wald/hausman type otherwise) -- deviance/            *
	 * deviance_below need no test() at all, vce(conventional) alone       *
	 * already posts them.                                                 */
	local basespec estimator(`estimator0') polynomial(`k0') `bwopt' ///
		limits(`L0' `Hlim0') `inv0' `mass0' `lg0' `nn0' allownegative ///
		deltamax(`deltamax') `rankopt' vce(conventional) nobias `spectestopt' nodrop

	/* every refit's window is inrange(zvar,lo,hi) from e(zlo)/e(zhi) --     *
	 * at the true cutoff that's exactly inrange(zvar,zlo0,zhi0), which      *
	 * this refit ALSO doubles as the bin-for-bin verification against       *
	 * Hstored (empirically confirmed to reproduce it exactly, together      *
	 * with nodrop above). If e(zlo)/e(zhi) are missing (an estimate that    *
	 * predates this feature), fall back to the literal if/in.               */
	local windowif0 "if inrange(`zvar', `zlo0', `zhi0')"
	if !`recenter' local windowif0 "`zif' `zin'"

	/* ------------------------------------------------------------------ *
	 *  Observed statistic, in one refit at the true cutoff.  By default  *
	 *  this doubles as the same-sample verification refit; clean skips   *
	 *  that check (the whole point of clean is that the loaded data is   *
	 *  deliberately NOT the data the model was fitted on).               *
	 * ------------------------------------------------------------------ */
	local rc0 = 0
	quietly {
		capture polbunch `zvarlist' `windowif0', cutoff(`cut0') `basespec'
		local rc0 = _rc
	}
	if `rc0' {
		di as error "could not re-fit the original specification on the currently loaded data"
		di as error "(rc=`rc0')."
		if "`clean'" == "" {
			di as error "Check that the loaded dataset is the one `=cond("`name'"=="","currently active","stored as `name'")' was fitted on."
		}
		exit 498
	}
	capture confirm matrix e(bins)
	if _rc {
		di as error "the true-cutoff refit did not produce e(bins)."
		exit 498
	}

	if "`clean'" == "" {
		tempname Hnow
		matrix `Hnow' = e(bins)
		local mismatch = (rowsof(`Hnow') != rowsof(`Hstored'))
		if !`mismatch' {
			local nr = rowsof(`Hstored')
			local nc = colsof(`Hstored')
			forvalues i = 1/`nr' {
				forvalues j = 1/`nc' {
					if abs(`Hnow'[`i',`j'] - `Hstored'[`i',`j']) > 1e-6 local mismatch = 1
				}
			}
		}
		if `mismatch' {
			di as error "the currently loaded dataset does not reproduce the stored model's"
			di as error "histogram bin-for-bin at cutoff(`cut0') bw(`bw0') -- is this the same data"
			di as error `"`=cond("`name'"=="","that produced the active estimate","`name' was fitted on")'?"'
			di as error "If you intend to test this specification against a genuinely different"
			di as error "(presumed kink-free) sample, use the clean option."
			exit 498
		}
	}

	capture local observed = _b[`target']
	if _rc {
		di as error "target(`target') is not a valid coefficient name."
		di as error "Choose a column of e(b), e.g. bunching:excess_mass (default),"
		di as error "bunching:number_bunchers, bunching:shift, bunching:marginal_response,"
		di as error "bunching:delta."
		exit 198
	}
	if missing(`observed') {
		di as error "target(`target') is missing/withheld at the true cutoff (weak identification"
		di as error "or an estimator that does not report it)."
		exit 498
	}
	capture local se_observed = _se[`target']
	if _rc | missing(`se_observed') | `se_observed' <= 0 {
		di as error "could not obtain a usable standard error for target(`target') at the true"
		di as error "cutoff (needed to studentize the permutation statistic)."
		exit 498
	}
	local observed_t = `observed' / `se_observed'

	/* spectests: pull the observed chi2/deviance values from this same   *
	 * true-cutoff refit -- e() still holds it at this point.             */
	local obs_omni  = .
	local obs_wald  = .
	local obs_haus  = .
	local obs_dev   = .
	local obs_devb  = .
	if "`spectests'" != "" {
		if `hadomni' local obs_omni = e(chi2_omnibus)
		if `hadwald' local obs_wald = e(chi2_wald)
		if `hadhaus' local obs_haus = e(chi2_hausman)
		if `haddev'  local obs_dev  = e(deviance)
		if `haddevb' local obs_devb = e(deviance_below)
	}

	/* ------------------------------------------------------------------ *
	 *  Candidate placebo-cutoff grid: step through the data range in     *
	 *  increments of step() (default the model's own bw), both          *
	 *  directions from the true cutoff, capped by maxcutoffs().          *
	 * ------------------------------------------------------------------ */
	local step = cond(`step'==0, `bw0', `step')

	/* every placebo brings its own recentered window (see above), so the   *
	 * candidate grid is no longer bounded by the target's literal if/in -- *
	 * it should span the full support of zvar in the loaded data. Falls    *
	 * back to the literal if/in only if e(zlo)/e(zhi) are unavailable      *
	 * (an estimate that predates this feature).                            */
	if `recenter' quietly summarize `zvar' `candif' `candin', meanonly
	else {
		local sumcond ""
		if `"`zifexpr'"' != "" & `"`candifexpr'"' != "" local sumcond "if (`zifexpr') & (`candifexpr')"
		else if `"`zifexpr'"' != "" local sumcond "if `zifexpr'"
		else if `"`candifexpr'"' != "" local sumcond "if `candifexpr'"
		quietly summarize `zvar' `sumcond' `zin' `candin', meanonly
	}
	local zdlo = r(min)
	local zdhi = r(max)

	/* Auto-detect candif/candin-induced data gaps and fold them into the     *
	 * SAME avoid()-style pre-filter as any explicit avoid() point, instead   *
	 * of only catching them reactively after a failed fit (see the candif   *
	 * classification in the fit loop below). Built on the RAW sorted data,  *
	 * not a bw0-phased bin grid: a "gap spell" is a maximal run of           *
	 * consecutive (in sorted z order) observations that all fail the        *
	 * caller's if/in, and its edges are the midpoints to the nearest         *
	 * surviving observations on each side -- exact and phase-independent,   *
	 * unlike binning at bw0, which would under-detect by up to a full bin    *
	 * width whenever the caller's excluded region (e.g. another kink) isn't  *
	 * phased the same as this model's own bin grid (a real case: two kinks   *
	 * in the same distribution need not sit on the same bw-spaced phase).    *
	 * A gap spell where NO observations ever existed on one side (grid      *
	 * edge) just uses its own boundary. Extra (lo,hi) pairs are appended     *
	 * after any explicit avoid() pairs, so a candidate whose window merely   *
	 * touches the hole is excluded before ever being fit -- exactly          *
	 * mirroring avoid()'s behaviour, just derived from the data instead of   *
	 * typed by the caller. Purely a data scan (preserve/restore, one sort,   *
	 * no model fits) -- a forvalues loop over SPELLS, not observations, so   *
	 * it stays cheap regardless of sample size. Only meaningful when a       *
	 * caller if/in was actually given, and only in the (common) recentered-  *
	 * window mode.                                                           */
	local navoidauto = 0
	if `recenter' & (`"`candifexpr'"' != "" | `"`candin'"' != "") {
		preserve
		quietly {
			keep if inrange(`zvar', `zdlo', `zdhi')
			tempvar ghascf gspell gzlo gzhi ghcnt
			gen byte `ghascf' = 0
			local gcfif ""
			if `"`candifexpr'"' != "" local gcfif "if (`candifexpr')"
			replace `ghascf' = 1 `gcfif' `candin'
			sort `zvar'
			gen long `gspell' = sum(`ghascf' != `ghascf'[_n-1])
			collapse (min) `gzlo'=`zvar' (max) `gzhi'=`zvar' (max) `ghcnt'=`ghascf', by(`gspell')
			sort `gspell'
			local ngs = _N
			forvalues i = 1/`ngs' {
				if `ghcnt'[`i'] == 0 {
					local glo = `gzlo'[`i']
					local ghi = `gzhi'[`i']
					if `i' > 1 local glo = (`gzhi'[`=`i'-1'] + `glo')/2
					if `i' < `ngs' local ghi = (`ghi' + `gzlo'[`=`i'+1'])/2
					local ++navoid
					local ++navoidauto
					local avoidlo`navoid' = `glo'
					local avoidhi`navoid' = `ghi'
				}
			}
		}
		restore
	}

	local maxperside = .
	if `maxcutoffs' != 0 local maxperside = floor(`maxcutoffs'/2)

	/* a candidate only belongs on the grid if its OWN full recentered      *
	 * window fits entirely within the available data -- not merely if     *
	 * the candidate cutoff itself does. Otherwise a candidate near the     *
	 * data's edge would silently get a truncated window (narrower than    *
	 * the target's own), rather than being excluded from consideration.   */
	local rmargin = cond(`recenter', `rightdist0', 0)
	local lmargin = cond(`recenter', `leftdist0',  0)

	/* window-overlap is filtered HERE, at grid construction, not left for   *
	 * the fit loop to discover: a candidate whose full recentered window    *
	 * would swallow the true excluded region (a real excess-mass spike,     *
	 * not "no data" -- see CLAUDE.md) never belongs in the candidate set at  *
	 * all. maxperside caps the number of GENUINELY VIABLE candidates per     *
	 * side, not the raw step count -- k keeps advancing (for free; this is   *
	 * pure arithmetic, no fit) past overlapping steps, on the chance a       *
	 * later one clears the excluded region and becomes viable again.        *
	 * Same right-closed/left-open asymmetry as elsewhere: lo_k>hi0 must be   *
	 * strict (hi0 itself belongs to the excluded bin).                      */
	local upcands ""
	local nup = 0
	local noverlap = 0
	local navoidoverlap = 0
	local k = 1
	while (`cut0' + `k'*`step' + `rmargin' <= `zdhi') & (`maxperside'>=. | `nup'<`maxperside') {
		local cand_k = `cut0' + `k'*`step'
		local iswoverlap = 0
		local isavoid = 0
		if `recenter' {
			local hi_k = `cand_k' + `rightdist0'
			local lo_k = `cand_k' - `leftdist0'
			local iswoverlap = !(`hi_k' <= `lo0' | `lo_k' > `hi0')
			if !`iswoverlap' & `navoid' > 0 {
				forvalues j = 1/`navoid' {
					if !(`hi_k' <= `avoidlo`j'' | `lo_k' > `avoidhi`j'') {
						local isavoid = 1
						continue, break
					}
				}
			}
		}
		if `iswoverlap' local ++noverlap
		else if `isavoid' local ++navoidoverlap
		else {
			local upcands "`upcands' `cand_k'"
			local ++nup
		}
		local ++k
	}
	local downcands ""
	local ndown = 0
	local k = 1
	while (`cut0' - `k'*`step' - `lmargin' >= `zdlo') & (`maxperside'>=. | `ndown'<`maxperside') {
		local cand_k = `cut0' - `k'*`step'
		local iswoverlap = 0
		local isavoid = 0
		if `recenter' {
			local hi_k = `cand_k' + `rightdist0'
			local lo_k = `cand_k' - `leftdist0'
			local iswoverlap = !(`hi_k' <= `lo0' | `lo_k' > `hi0')
			if !`iswoverlap' & `navoid' > 0 {
				forvalues j = 1/`navoid' {
					if !(`hi_k' <= `avoidlo`j'' | `lo_k' > `avoidhi`j'') {
						local isavoid = 1
						continue, break
					}
				}
			}
		}
		if `iswoverlap' local ++noverlap
		else if `isavoid' local ++navoidoverlap
		else {
			local downcands "`downcands' `cand_k'"
			local ++ndown
		}
		local ++k
	}
	local ncand = `nup' + `ndown'
	if `ncand' == 0 {
		di as error "no placebo cutoffs fit within the data range at step(`step')"
		di as error "without overlapping the true excluded region (`noverlap' candidate(s)) or an"
		di as error "avoid() region (`navoidoverlap' candidate(s))."
		exit 498
	}

	tempname draws
	matrix `draws' = J(`ncand', 4, .)
	tempname stdraws
	if "`spectests'" != "" matrix `stdraws' = J(`ncand', 7, .)
	/* stdraws has no status column of its own -- it shares placebo_draws'   *
	 * (used/failed), since a candidate invalid for excess_mass is equally   *
	 * invalid for its spec-test statistics. Overlap no longer appears as a  *
	 * status value here -- overlapping candidates were already excluded     *
	 * above and never enter this matrix at all. A candidate that converges  *
	 * but fails only ITS test computation just leaves the specific          *
	 * chi2/deviance column missing, filtered out at aggregation. Column 7   *
	 * (only meaningful when a chi2 test dropped norankcheck, nchi2types>0)  *
	 * tracks e(polynomial) per placebo, since that's the one thing          *
	 * norankcheck was otherwise guaranteeing stays fixed.                   */
	local i = 0
	foreach cand_k of local upcands {
		local ++i
		matrix `draws'[`i',1] = `cand_k'
		if "`spectests'" != "" matrix `stdraws'[`i',1] = `cand_k'
	}
	foreach cand_k of local downcands {
		local ++i
		matrix `draws'[`i',1] = `cand_k'
		if "`spectests'" != "" matrix `stdraws'[`i',1] = `cand_k'
	}

	/* ------------------------------------------------------------------ *
	 *  Placebo loop -- polbunch itself preserve/restores the loaded      *
	 *  data on every call (including on error), so this is safe to run   *
	 *  directly on the live dataset with no frame/preserve wrapper here. *
	 * ------------------------------------------------------------------ */
	local nused = 0
	local nfailed = 0
	local ncandif = 0
	/* noverlap NOT reset here -- it already holds the grid-construction    *
	 * count (candidates excluded before ever reaching this loop); the      *
	 * fit loop below only ADDS to it, for the window_recentered==0         *
	 * fallback where overlap can only be discovered after fitting.        */
	if "`nodots'" == "" nois _dots 0, title("Refitting at placebo cutoffs...") reps(`ncand')
	forvalues i = 1/`ncand' {
		local cand = `draws'[`i',1]
		local rc_i = 0
		local val = .
		local se_i = .
		local lok = .
		local hik = .
		local pl_omni = .
		local pl_wald = .
		local pl_haus = .
		local pl_dev  = .
		local pl_devb = .
		local pl_poly = .
		local windowif_i "`zif' `zin'"
		local lo_i = .
		local hi_i = .
		if `recenter' {
			local lo_i = `cand' - `leftdist0'
			local hi_i = `cand' + `rightdist0'
			local wcond_i "inrange(`zvar', `lo_i', `hi_i')"
			if `"`candifexpr'"' != "" local wcond_i "`wcond_i' & (`candifexpr')"
			local windowif_i "if `wcond_i' `candin'"
		}
		else if `"`candifexpr'"' != "" | `"`candin'"' != "" {
			local wcond_i "`zifexpr'"
			if `"`candifexpr'"' != "" {
				if `"`wcond_i'"' != "" local wcond_i "(`wcond_i') & (`candifexpr')"
				else local wcond_i "`candifexpr'"
			}
			local windowif_i = cond(`"`wcond_i'"'=="", "", "if `wcond_i'")
			local windowif_i "`windowif_i' `zin' `candin'"
		}
		/* window-level overlap, checked BEFORE fitting: even when this      *
		 * candidate's own tiny excluded donut is nowhere near the true      *
		 * excluded region, its wider REFERENCE window can still reach into  *
		 * it -- and the true excluded region is not just "no data" (there   *
		 * is no missing-mass hole under the iso-elastic model, see          *
		 * CLAUDE.md), it is a genuine excess-mass SPIKE, a real structural  *
		 * discontinuity. Folding that spike into another cutoff's ordinary  *
		 * reference-region fit as if it were unremarkable data would        *
		 * corrupt that placebo's counterfactual for reasons having nothing  *
		 * to do with whether there is a real effect AT the placebo's own    *
		 * location. Skips the fit entirely -- cheaper, and the widened      *
		 * grid (full data support, not just the target's literal window)    *
		 * makes far more candidates reach this than before.                 *
		 *                                                                    *
		 * Asymmetric strictness, NOT a typo: polbunch bins are right-closed, *
		 * left-open (z==cutoff belongs to the bin BELOW it), so the true     *
		 * excluded region (lo0,hi0] itself EXCLUDES lo0 but INCLUDES hi0.    *
		 * hi_i<=lo0 (weak) correctly says "no overlap" when the window ends  *
		 * exactly at lo0, since lo0 belongs to the bin below, not the        *
		 * excluded one. But lo_i==hi0 is NOT the same kind of touch: hi0 IS  *
		 * part of the excluded bin -- exactly where a real bunching spike    *
		 * concentrates as a discrete atom, not an infinitesimal continuous   *
		 * density, so this boundary can carry real, substantial mass. A      *
		 * window starting AT hi0 must count as overlapping, hence lo_i>hi0   *
		 * (strict), not lo_i>=hi0.                                           */
		local windowoverlap = 0
		if `recenter' {
			local windowoverlap = !(`hi_i' <= `lo0' | `lo_i' > `hi0')
			if !`windowoverlap' & `navoid' > 0 {
				forvalues j = 1/`navoid' {
					if !(`hi_i' <= `avoidlo`j'' | `lo_i' > `avoidhi`j'') {
						local windowoverlap = 1
						continue, break
					}
				}
			}
		}
		if `windowoverlap' {
			matrix `draws'[`i',4] = -1
			local ++noverlap
		}
		else {
		quietly {
			capture polbunch `zvarlist' `windowif_i', cutoff(`cand') `basespec'
			local rc_i = _rc
			if !`rc_i' {
				capture local val = _b[`target']
				capture local se_i = _se[`target']
				local lok = e(lower_limit)
				local hik = e(upper_limit)
				if missing(`lok') | missing(`hik') {
					local lok = e(zL_excl_orig)
					local hik = e(zH_excl_orig)
				}
				if "`spectests'" != "" {
					if `hadomni' local pl_omni = e(chi2_omnibus)
					if `hadwald' local pl_wald = e(chi2_wald)
					if `hadhaus' local pl_haus = e(chi2_hausman)
					if `haddev'  local pl_dev  = e(deviance)
					if `haddevb' local pl_devb = e(deviance_below)
					if `nchi2types' > 0 local pl_poly = e(polynomial)
				}
			}
		}
		if `rc_i' | missing(`val') | missing(`se_i') | `se_i'<=0 | missing(`lok') | missing(`hik') {
			/* Before counting this as a genuine non-convergence: if the      *
			 * caller passed its own if/in (candifexpr/candin) to this        *
			 * command -- e.g. to carve out ANOTHER real kink's excluded      *
			 * region -- a candidate whose window straddles that hole has no  *
			 * data to identify one side of the fit, and fails for a reason   *
			 * that has nothing to do with whether there's bunching at this   *
			 * placebo. Re-run the identical window WITHOUT the caller's      *
			 * if/in: if that succeeds, the original failure was caused       *
			 * purely by the caller's own intentional restriction, not by     *
			 * this being a bad placebo location -- drop it from the grid     *
			 * entirely (see the Mata filter below) rather than reporting it  *
			 * as a failed draw. This is a classification probe only; its     *
			 * fitted value is discarded either way, since it would still be  *
			 * contaminated by whatever the caller's if/in was excluding.     */
			local iscandif = 0
			if `"`candifexpr'"' != "" | `"`candin'"' != "" {
				local rc_i2  = 1
				local val2   = .
				local se_i2  = .
				local lok2   = .
				local hik2   = .
				local windowif_i2 "`zif' `zin'"
				if `recenter' local windowif_i2 "if inrange(`zvar', `lo_i', `hi_i')"
				else if `"`zifexpr'"' != "" local windowif_i2 "if `zifexpr' `zin'"
				quietly {
					capture polbunch `zvarlist' `windowif_i2', cutoff(`cand') `basespec'
					local rc_i2 = _rc
					if !`rc_i2' {
						capture local val2 = _b[`target']
						capture local se_i2 = _se[`target']
						local lok2 = e(lower_limit)
						local hik2 = e(upper_limit)
						if missing(`lok2') | missing(`hik2') {
							local lok2 = e(zL_excl_orig)
							local hik2 = e(zH_excl_orig)
						}
					}
				}
				local iscandif = !`rc_i2' & !missing(`val2') & !missing(`se_i2') & `se_i2'>0 & !missing(`lok2') & !missing(`hik2')
			}
			if `iscandif' {
				matrix `draws'[`i',4] = -2
				local ++ncandif
			}
			else {
				matrix `draws'[`i',4] = 0
				local ++nfailed
			}
		}
		else {
			/* same right-closed/left-open asymmetry as the window-overlap    *
			 * check above: hik<=lo0 (weak) is correct, lok>hi0 must be       *
			 * strict, since hi0 itself belongs to the excluded bin.          */
			local overlaps = !(`hik' <= `lo0' | `lok' > `hi0')
			if !`overlaps' & `navoid' > 0 {
				forvalues j = 1/`navoid' {
					if !(`hik' <= `avoidlo`j'' | `lok' > `avoidhi`j'') {
						local overlaps = 1
						continue, break
					}
				}
			}
			if `overlaps' {
				matrix `draws'[`i',4] = -1
				local ++noverlap
			}
			else {
				matrix `draws'[`i',2] = `val'
				matrix `draws'[`i',3] = `se_i'
				matrix `draws'[`i',4] = 1
				local ++nused
				if "`spectests'" != "" {
					matrix `stdraws'[`i',2] = `pl_omni'
					matrix `stdraws'[`i',3] = `pl_wald'
					matrix `stdraws'[`i',4] = `pl_haus'
					matrix `stdraws'[`i',5] = `pl_dev'
					matrix `stdraws'[`i',6] = `pl_devb'
					matrix `stdraws'[`i',7] = `pl_poly'
				}
			}
		}
		}
		if "`nodots'" == "" nois _dots `i' 0
	}

	/* Drop candidates whose only failure was the caller's own if/in         *
	 * (status -2, set above) -- these never belonged on the grid, they      *
	 * just couldn't be recognised as such until after an attempted fit.     *
	 * Removed from placebo_draws (and spectest_draws, same row selection)   *
	 * entirely, not merely relabelled, so they don't inflate "attempted"/   *
	 * "failed" counts in the report below; ncand is refreshed to match.     */
	if `ncandif' > 0 {
		mata: _pbx_keep = st_matrix("`draws'")[.,4] :!= -2
		mata: st_matrix("`draws'", select(st_matrix("`draws'"), _pbx_keep))
		if "`spectests'" != "" mata: st_matrix("`stdraws'", select(st_matrix("`stdraws'"), _pbx_keep))
		mata: mata drop _pbx_keep
		local ncand = rowsof(`draws')
	}

	/* sort ascending by cutoff (column 1) -- rows were filled up-candidates  *
	 * first (increasing), then down-candidates (decreasing), so the raw     *
	 * fill order is not monotonic. Apply the SAME row permutation to        *
	 * stdraws so the two matrices stay aligned candidate-for-candidate.     */
	mata: st_matrix("_pbx_sortorder", order(st_matrix("`draws'"), 1))
	mata: st_matrix("`draws'", st_matrix("`draws'")[st_matrix("_pbx_sortorder"), .])
	if "`spectests'" != "" mata: st_matrix("`stdraws'", st_matrix("`stdraws'")[st_matrix("_pbx_sortorder"), .])

	matrix colnames `draws' = cutoff value se status
	if "`spectests'" != "" matrix colnames `stdraws' = cutoff chi2_omnibus chi2_wald chi2_hausman deviance deviance_below polynomial

	/* ------------------------------------------------------------------ *
	 *  Aggregate: exact permutation p-value ("+1" correction) on the     *
	 *  STUDENTIZED statistic (value/se), plus descriptive stats of both  *
	 *  the raw and studentized placebo distributions.                    *
	 * ------------------------------------------------------------------ */
	local pval = .
	local vmean = .
	local vsd   = .
	local vmin  = .
	local vmax  = .
	local vlb   = .
	local vub   = .
	local tmean = .
	local tsd   = .
	local tlb   = .
	local tub   = .
	if `nused' > 0 {
		local cnt = 0
		forvalues i = 1/`ncand' {
			if `draws'[`i',4] == 1 {
				local t = `draws'[`i',2] / `draws'[`i',3]
				if "`onesided'" == "" {
					if abs(`t') >= abs(`observed_t') local ++cnt
				}
				else {
					if `t' >= `observed_t' local ++cnt
				}
			}
		}
		local pval = (1 + `cnt') / (1 + `nused')

		tempname fr2
		frame create `fr2'
		frame `fr2' {
			quietly {
				svmat double `draws', names(_pd)
				gen double _t = _pd2/_pd3 if _pd4==1
				summarize _pd2 if _pd4==1
				local vmean = r(mean)
				local vsd   = r(sd)
				local vmin  = r(min)
				local vmax  = r(max)
				summarize _t if _pd4==1
				local tmean = r(mean)
				local tsd   = r(sd)
				local a2 = (100 - `level')/2
				_pctile _pd2 if _pd4==1, p(`a2' `=100-`a2'')
				local vlb = r(r1)
				local vub = r(r2)
				_pctile _t if _pd4==1, p(`a2' `=100-`a2'')
				local tlb = r(r1)
				local tub = r(r2)
			}
		}
		frame drop `fr2'
	}

	/* ------------------------------------------------------------------ *
	 *  spectests aggregation: one-sided permutation p-value per detected *
	 *  test (chi2/deviance are non-negative, larger = more extreme --    *
	 *  the "onesided" option governs the excess_mass p-value only, not   *
	 *  these). Uses the same placebo_draws status (col 4) for the        *
	 *  used/overlap/failed classification; a candidate that converged    *
	 *  but whose OWN test computation came back missing is dropped only  *
	 *  from that specific test's count, not from the others.             *
	 * ------------------------------------------------------------------ */
	local p_omni_perm  = .
	local p_wald_perm  = .
	local p_hausman_perm = .
	local p_deviance_perm = .
	local p_deviance_below_perm = .
	local n_omni  = 0
	local n_wald  = 0
	local n_haus  = 0
	local n_dev   = 0
	local n_devb  = 0
	local degree_drift = 0
	if "`spectests'" != "" {
		local cnt_omni = 0
		local cnt_wald = 0
		local cnt_haus = 0
		local cnt_dev  = 0
		local cnt_devb = 0
		forvalues i = 1/`ncand' {
			if `draws'[`i',4] == 1 {
				if `hadomni' & !missing(`stdraws'[`i',2]) {
					local ++n_omni
					if `stdraws'[`i',2] >= `obs_omni' local ++cnt_omni
				}
				if `hadwald' & !missing(`stdraws'[`i',3]) {
					local ++n_wald
					if `stdraws'[`i',3] >= `obs_wald' local ++cnt_wald
				}
				if `hadhaus' & !missing(`stdraws'[`i',4]) {
					local ++n_haus
					if `stdraws'[`i',4] >= `obs_haus' local ++cnt_haus
				}
				if `haddev' & !missing(`stdraws'[`i',5]) {
					local ++n_dev
					if `stdraws'[`i',5] >= `obs_dev' local ++cnt_dev
				}
				if `haddevb' & !missing(`stdraws'[`i',6]) {
					local ++n_devb
					if `stdraws'[`i',6] >= `obs_devb' local ++cnt_devb
				}
			}
		}
		if `hadomni' & `n_omni' > 0 local p_omni_perm = (1 + `cnt_omni') / (1 + `n_omni')
		if `hadwald' & `n_wald' > 0 local p_wald_perm = (1 + `cnt_wald') / (1 + `n_wald')
		if `hadhaus' & `n_haus' > 0 local p_hausman_perm = (1 + `cnt_haus') / (1 + `n_haus')
		if `haddev'  & `n_dev'  > 0 local p_deviance_perm = (1 + `cnt_dev') / (1 + `n_dev')
		if `haddevb' & `n_devb' > 0 local p_deviance_below_perm = (1 + `cnt_devb') / (1 + `n_devb')

		/* degree drift: norankcheck was dropped to get the chi2 tests at  *
		 * all (see above), so a placebo's polynomial can in principle     *
		 * come back lower than the target's own k0 -- count how often.   */
		local degree_drift = 0
		if `nchi2types' > 0 {
			forvalues i = 1/`ncand' {
				if `draws'[`i',4] == 1 & !missing(`stdraws'[`i',7]) & `stdraws'[`i',7] != `k0' {
					local ++degree_drift
				}
			}
		}
	}

	/* ------------------------------------------------------------------ *
	 *  Report                                                             *
	 * ------------------------------------------------------------------ */
	local sidelab = cond("`onesided'"=="","two-sided","one-sided")
	di ""
	if "`clean'" == "" {
		di as text "Placebo-cutoff permutation test of {res:`tlab'}{txt} for a polbunch estimate"
	}
	else {
		di as text "Placebo-cutoff permutation test of {res:`tlab'}{txt} (CLEAN: cross-sample mode)"
		di as text "the currently loaded data was NOT verified to be the sample `=cond("`name'"=="","the active estimate","`name'")' was fitted on."
	}
	di as text "{hline 74}"
	di as text %-40s "true cutoff" _col(66) as res %12.0g `cut0'
	di as text %-40s "observed `tlab' (se)" _col(55) as res %10.4f `observed' as text " (" as res %6.4f `se_observed' as text ")"
	di as text %-40s "observed `tlab', studentized" _col(66) as res %10.4f `observed_t'
	di as text "{hline 74}"
	di as text %-40s "excluded before fitting (window would overlap true region)" _col(66) as res %10.0f `noverlap'
	if `navoid' > 0 {
		di as text %-40s "excluded before fitting (window overlaps avoid()/if-in gap)" _col(66) as res %10.0f `navoidoverlap'
	}
	if `ncandif' > 0 {
		di as text %-40s "excluded (window overlaps your own if/in restriction)" _col(66) as res %10.0f `ncandif'
	}
	di as text %-40s "placebo cutoffs attempted (remaining)" _col(66) as res %10.0f `ncand'
	di as text %-40s "  used (converged)" _col(66) as res %10.0f `nused'
	di as text %-40s "  failed to converge / no usable SE" _col(66) as res %10.0f `nfailed'
	if `nused' > 0 {
		di as text "{hline 74}"
		di as text %-40s "placebo `tlab' (raw): mean / sd" _col(45) as res %14.4f `vmean' as text " /" as res %13.4f `vsd'
		di as text %-40s "placebo `tlab' (raw): min / max" _col(45) as res %14.4f `vmin' as text " /" as res %13.4f `vmax'
		di as text %-40s "`level'% range of raw placebo draws" _col(45) as res %14.4f `vlb' as text " /" as res %13.4f `vub'
		di as text %-40s "placebo `tlab' (studentized): mean / sd" _col(45) as res %14.4f `tmean' as text " /" as res %13.4f `tsd'
		di as text %-40s "`level'% range of studentized draws" _col(45) as res %14.4f `tlb' as text " /" as res %13.4f `tub'
		di as text "{hline 74}"
		di as text %-40s "H0: no true bunching (`sidelab' permutation p)" _col(66) as res %10.4f `pval'
	}
	else {
		di as text "No placebo cutoff produced a usable draw -- p-value not computed."
	}
	if "`spectests'" != "" & (`hadomni' | `hadwald' | `hadhaus' | `haddev' | `haddevb') {
		di as text "{hline 74}"
		di as text "spectests: permutation p-values for the target's own tests (one-sided," ///
			" larger = more extreme)"
		if `hadomni'  di as text %-40s "Omnibus chi2 = " %8.3f `obs_omni'  _col(66) as res %10.4f `p_omni_perm'
		if `hadwald'  di as text %-40s "Wald chi2 (off-label diagnostic) = " %8.3f `obs_wald'  _col(66) as res %10.4f `p_wald_perm'
		if `hadhaus'  di as text %-40s "Hausman chi2 = " %8.3f `obs_haus'  _col(66) as res %10.4f `p_hausman_perm'
		if `haddev'   di as text %-40s "Deviance = " %8.3f `obs_dev'   _col(66) as res %10.4f `p_deviance_perm'
		if `haddevb'  di as text %-40s "Deviance, below cutoff only = " %8.3f `obs_devb'  _col(66) as res %10.4f `p_deviance_below_perm'
		if `nchi2types' > 0 {
			di as text "note: a chi2 test was requested, so norankcheck was dropped (polbunch forces" ///
				" test(none) under norankcheck) -- placebo polynomial degree differs from the"
			di as text "target's own (`k0') in " as res `degree_drift' as text " of `nused' used placebo cutoffs."
		}
	}
	di as text "{hline 74}"
	di as text "step(`step'), deterministic grid (no resampling, no seed); permutation test is on"
	di as text "the studentized statistic (target / analytic SE), not raw `tlab'."
	if `recenter' {
		di as text "note: the target's if/in was purely window-defining, so every placebo's own"
		di as text "estimation window was RE-CENTERED on its own cutoff (halfwidths " ///
			as res %6.4f `leftdist0' as text " / " as res %6.4f `rightdist0' as text " left/right)"
		di as text "rather than reusing the target's literal numeric window bounds."
	}
	if `"`candif'"' != "" | `"`candin'"' != "" {
		di as text "note: candidate data further restricted by " as res "`candif' `candin'" ///
			as text " (candidate grid and placebo windows only, not the true-cutoff refit)."
	}
	if `navoiduser' > 0 {
		di as text "note: avoid() kept " as res `navoiduser' as text " extra location(s) out of the candidate" ///
			" grid, using the target's own excluded-region half-widths recentred at each."
	}
	if `navoidauto' > 0 {
		di as text "note: your if/in emptied " as res `navoidauto' as text " contiguous data gap(s) -- candidates" ///
			" whose window overlaps any of them were excluded before fitting, same as avoid()."
	}

	/* ------------------------------------------------------------------ *
	 *  Return                                                             *
	 * ------------------------------------------------------------------ */
	if "`curkind'" == "hold" _estimates unhold `_cur'
	else {
		quietly estimates restore `_cur'
		capture estimates drop `_cur'
	}

	return matrix placebo_draws = `draws'
	if `nused' > 0 {
		return scalar p              = `pval'
		return scalar placebo_mean   = `vmean'
		return scalar placebo_sd     = `vsd'
		return scalar placebo_min    = `vmin'
		return scalar placebo_max    = `vmax'
		return scalar placebo_lb     = `vlb'
		return scalar placebo_ub     = `vub'
		return scalar placebo_t_mean = `tmean'
		return scalar placebo_t_sd   = `tsd'
		return scalar placebo_t_lb   = `tlb'
		return scalar placebo_t_ub   = `tub'
	}
	return scalar observed        = `observed'
	return scalar observed_se     = `se_observed'
	return scalar observed_t      = `observed_t'
	return scalar n_placebo       = `nused'
	return scalar n_placebo_tried = `ncand'
	return scalar cutoffs_skipped = `nfailed'
	return scalar cutoffs_overlap = `noverlap'
	return scalar cutoffs_avoid   = `navoidoverlap'
	return scalar avoid_user      = `navoiduser'
	return scalar avoid_auto      = `navoidauto'
	return scalar cutoffs_candif  = `ncandif'
	return scalar step            = `step'
	return scalar bw              = `bw0'
	return scalar cutoff_true     = `cut0'
	return scalar onesided        = ("`onesided'" != "")
	return scalar clean           = ("`clean'" != "")
	return scalar level           = `level'
	return scalar spectests       = ("`spectests'" != "")
	return scalar window_recentered = `recenter'
	return scalar window_leftdist   = `leftdist0'
	return scalar window_rightdist  = `rightdist0'
	if "`spectests'" != "" {
		return matrix spectest_draws = `stdraws'
		return scalar spectest_degree_drift = `degree_drift'
		if `hadomni' {
			return scalar observed_omnibus = `obs_omni'
			return scalar p_omnibus_perm   = `p_omni_perm'
		}
		if `hadwald' {
			return scalar observed_wald = `obs_wald'
			return scalar p_wald_perm   = `p_wald_perm'
		}
		if `hadhaus' {
			return scalar observed_hausman = `obs_haus'
			return scalar p_hausman_perm   = `p_hausman_perm'
		}
		if `haddev' {
			return scalar observed_deviance = `obs_dev'
			return scalar p_deviance_perm   = `p_deviance_perm'
		}
		if `haddevb' {
			return scalar observed_deviance_below = `obs_devb'
			return scalar p_deviance_below_perm   = `p_deviance_below_perm'
		}
	}
	return local target "`target'"
	return local name "`name'"
end
