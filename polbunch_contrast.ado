*! polbunch_contrast version date 20260917
* Author: Martin Eckhoff Andresen
* This program is part of the polbunch package.
*
* Postestimation: a paired binned-bootstrap test of the difference in a
* parameter (default the implied elasticity) between TWO stored polbunch
* estimates that were fitted on the SAME histogram.  Both specifications are
* re-fitted on every resampled histogram, so the bootstrap distribution of the
* difference carries the dependence of the two estimators (they share the bin
* counts).  Nothing needs to be in memory -- the histogram is read from e(bins)
* of the stored estimates.
*
*   polbunch_contrast nameA nameB [, target(coef) reps(#) seed(#) ///
*                                   boottype(multinomial) level(#) nodots ]
*
* Returns r(diff), r(se), r(z), r(p), r(p_bootstrap), r(corr), r(ci_ll),
* r(ci_ul), r(value_A), r(value_B), r(reps), r(reps_used), r(target),
* r(spec_A), r(spec_B), r(name_A), r(name_B), r(warning).

cap program drop polbunch_contrast
program define polbunch_contrast, rclass
	version 16.0

	syntax namelist(min=2 max=2) [, ///
		TARGET(string) ///
		reps(integer 999) ///
		SEED(string) ///
		BOOTtype(string) ///
		Level(cilevel) ///
		NODOTS ]

	gettoken nameA nameB : namelist
	local nameA = strtrim("`nameA'")
	local nameB = strtrim("`nameB'")

	/* which coefficient to contrast -- any element of e(b), e.g.
	   bunching:shift, bunching:number_bunchers, bunching:delta.  Default the
	   implied elasticity. */
	if "`target'" == "" local target "bunching:elasticity"
	local target = strtrim("`target'")
	local tlab : subinstr local target "bunching:" "", all
	local tlab : subinstr local tlab "_" " ", all

	if `reps' < 2 {
		di as error "reps() must be at least 2."
		exit 198
	}
	if "`boottype'" == "" local boottype multinomial
	local boottype = strlower("`boottype'")
	if !inlist("`boottype'","multinomial") {
		di as error "boottype() must be multinomial.  residual / wild resample around"
		di as error "a per-model fitted mean and are not supported here; use the contrast option"
		di as error "inside polbunch for those."
		exit 198
	}
	if "`seed'" != "" {
		capture set seed `seed'
		if _rc {
			di as error "invalid seed()."
			exit 198
		}
	}

	/* hold whatever estimates are active now; restored on exit / error */
	tempname _cur
	_estimates hold `_cur', restore nullok

	/* ------------------------------------------------------------------ *
	 *  Read the two stored estimates                                      *
	 * ------------------------------------------------------------------ */
	tempname HA HB
	local j 0
	foreach sname in `nameA' `nameB' {
		local ++j
		local tag = cond(`j'==1,"A","B")

		capture estimates restore `sname'
		if _rc {
			di as error "could not restore `sname' -- use {cmd:estimates store `sname'} after fitting."
			exit 198
		}
		if "`e(cmd)'" != "polbunch" {
			di as error "`sname' is not a polbunch estimate."
			exit 198
		}
		if !inlist(e(estimator),0,1,2,3,4) {
			di as error "`sname' has an unrecognised e(estimator)."
			exit 198
		}
		capture confirm matrix e(bins)
		if _rc {
			di as error "`sname' has no e(bins); it was fitted with an older polbunch."
			di as error "Re-fit `sname' and {cmd:estimates store} it again."
			exit 198
		}
		tempname _el`tag'
		capture scalar `_el`tag'' = _b[`target']
		if _rc {
			di as error "target(`target') is not in `sname''s coefficient vector."
			di as error "Choose a column of e(b), e.g. bunching:elasticity, bunching:shift,"
			di as error "bunching:number_bunchers, bunching:excess_mass, bunching:delta."
			exit 198
		}
		if missing(`_el`tag'') {
			di as error "`sname' has `target' withheld / not identified (weak delta or no tax rates)."
			exit 198
		}

		matrix `H`tag'' = e(bins)

		local est`tag'    = e(estimator)
		local k`tag'      = e(polynomial)
		local cut`tag'    = e(cutoff_orig)
		local bw`tag'     = e(bw)
		local lo`tag'     = e(lower_limit)
		local hi`tag'     = e(upper_limit)
		local con`tag'    = e(constant)
		local nsp`tag'    = e(nosplit)
		local t0`tag'     = e(t0)
		local t1`tag'     = e(t1)
		local islog`tag'  = (e(log) == 1)
		local nonorm`tag' = ("`e(normalize)'" == "nonormalize")
	}

	/* ------------------------------------------------------------------ *
	 *  The two histograms must be the same data, or one must be a       *
	 *  bin-for-bin subset of the other (e.g. B fit on a narrower        *
	 *  z-range of the same underlying data, same bin grid).  Either     *
	 *  way, the resample below draws from whichever histogram is the    *
	 *  superset, and each spec is refit on only the bins it actually    *
	 *  used.                                                             *
	 * ------------------------------------------------------------------ */
	local nA = rowsof(`HA')
	local nB = rowsof(`HB')
	if abs(`bwA' - `bwB') > 1e-8*max(abs(`bwA'),abs(`bwB')) {
		di as error "the two estimates use different bin widths (`bwA' vs `bwB') -- their"
		di as error "histograms cannot be aligned; polbunch_contrast requires the same"
		di as error "underlying bin grid."
		exit 198
	}

	tempname Hbig Hsmall
	if `nA' >= `nB' {
		matrix `Hbig'   = `HA'
		matrix `Hsmall' = `HB'
		local bigtag   A
		local smalltag B
	}
	else {
		matrix `Hbig'   = `HB'
		matrix `Hsmall' = `HA'
		local bigtag   B
		local smalltag A
	}
	local nbig   = rowsof(`Hbig')
	local nsmall = rowsof(`Hsmall')

	/* locate the smaller histogram inside the bigger one by matching bin
	   midpoints; every row of the smaller must appear, in order, as a
	   contiguous run in the bigger, with matching counts, for the pair
	   to qualify as "same data" (nsmall==nbig) or a nested subset. */
	local nested 1
	local offset .
	forvalues i = 1/`nsmall' {
		local found .
		forvalues j = 1/`nbig' {
			if abs(`Hbig'[`j',2] - `Hsmall'[`i',2]) < `bwA'*1e-4 {
				local found `j'
				continue, break
			}
		}
		if `found' >= . {
			local nested 0
			continue, break
		}
		if reldif(`Hbig'[`found',1], `Hsmall'[`i',1]) > 1e-7 & ///
			abs(`Hbig'[`found',1] - `Hsmall'[`i',1]) > 1e-6 {
			local nested 0
			continue, break
		}
		if `i' == 1 local offset = `found'
		else if `found' != `offset' + `i' - 1 {
			local nested 0
			continue, break
		}
	}
	if !`nested' {
		di as error "the two estimates were fitted on different histograms (counts or bin"
		di as error "positions differ, and neither is a bin-for-bin subset of the other) --"
		di as error "polbunch_contrast requires the same data, or one estimate to be fit on"
		di as error "a narrower range of the other's data."
		exit 198
	}
	local samerange = (`nsmall' == `nbig')

	/* if the specs use different ranges, the narrower one's polbunch     *
	 * refit must see only its own bins (they are rows `offset'..         *
	 * `offset'+`nsmall'-1 of the resampled `Hbig' frame built below).    */
	local restrictA ""
	local restrictB ""
	if !`samerange' {
		local _restrict "if inrange(_n,`offset',`=`offset'+`nsmall'-1')"
		if "`smalltag'" == "A" local restrictA "`_restrict'"
		else                   local restrictB "`_restrict'"
	}

	/* ------------------------------------------------------------------ *
	 *  Recover limits(L H) from the stored excluded-region edges.        *
	 *  polbunch has two conventions: if the cutoff falls on a bin EDGE,  *
	 *  limits(L H) excludes L whole bins below and H above with no       *
	 *  cutoff bin (lower_limit = cutoff - L*bw);  if the cutoff falls    *
	 *  INSIDE a bin, that bin plus L below and H above are excluded      *
	 *  (lower_limit = cutoff - (L + 0.5)*bw).  Detect which from the     *
	 *  offset of the cutoff within the bin grid.                         *
	 * ------------------------------------------------------------------ */
	foreach tag in A B {
		local _cut = cond("`tag'"=="A",`cutA',`cutB')
		local _bw  = cond("`tag'"=="A",`bwA',`bwB')
		local _lo  = cond("`tag'"=="A",`loA',`loB')
		local _hi  = cond("`tag'"=="A",`hiA',`hiB')
		local _m0  = `H`tag''[1,2]
		local _frac = mod((`_cut' - `_m0')/`_bw', 1)
		if `_frac' > 0.5 local _frac = `_frac' - 1
		local _edge`tag' = (abs(abs(`_frac') - 0.5) < 0.25)
		if `_edge`tag'' {
			local L_`tag'   = round((`_cut' - `_lo')/`_bw')
			local Hlim_`tag' = round((`_hi' - `_cut')/`_bw')
		}
		else {
			local L_`tag'   = round((`_cut' - `_lo')/`_bw' - 0.5)
			local Hlim_`tag' = round((`_hi' - `_cut')/`_bw' - 0.5)
		}
	}
	local LA = `L_A'
	local HAlim = `Hlim_A'
	local LB = `L_B'
	local HBlim = `Hlim_B'

	local warn ""
	if `kA' != `kB' & !inlist(4,`estA',`estB') ///
	                                      local warn "`warn' polynomial order (`kA' vs `kB');"
	if `LA' != `LB' | `HAlim' != `HBlim'   local warn "`warn' excluded window (limits `LA' `HAlim' vs `LB' `HBlim');"
	if reldif(`cutA',`cutB') > 1e-8        local warn "`warn' cutoff (`cutA' vs `cutB');"
	if !`samerange'                        local warn "`warn' estimation range (`bigtag' uses `nbig' bins, `smalltag' a narrower `nsmall'-bin subset of the same data);"
	local warn = strtrim("`warn'")

	/* a pair that mixes different counterfactual structures -- the free
	   two-sided fit (0), the restricted model (3), or the Saez two-point
	   trapezoid (4) -- is contrasting estimands that need not share a plim
	   even under the model.  Flag it. */
	local mixnote ""
	if inlist(0,`estA',`estB') | inlist(4,`estA',`estB') ///
		local mixnote "the two estimators use different counterfactual models"

	/* ------------------------------------------------------------------ *
	 *  Canonical pre-binned polbunch spec strings                        *
	 * ------------------------------------------------------------------ */
	foreach tag in A B {
		local inv`tag'  = cond(`con`tag'',"constant","exact")
		local mass`tag' = cond(`nsp`tag'',"poolmass","splitmass")
		local lg`tag'   = cond(`islog`tag'',"log","")
		local nn`tag'   = cond(`nonorm`tag'',"nonormalize","")
		local H_`tag'   = `Hlim_`tag''
		/* norankcheck: hold the polynomial degree fixed at what the stored fit
		   used (e(polynomial) already reflects any automatic reduction), so the
		   resampled refits do not drift to a lower degree on awkward draws --
		   that would inflate the bootstrap variance relative to the stored
		   specification. */
		local spec`tag' estimator(`est`tag'') polynomial(`k`tag'') cutoff(`cut`tag'') ///
			t0(`t0`tag'') t1(`t1`tag'') limits(`L_`tag'' `H_`tag'') ///
			`inv`tag'' `mass`tag'' `lg`tag'' `nn`tag'' ///
			norankcheck vce(none) nobias test(none)

		if `est`tag'' == 4      local lab`tag' "est 4 (Saez trapezoid)"
		else if `est`tag'' == 0 local lab`tag' "est 0 (left data only), poly `k`tag''"
		else                    local lab`tag' "est `est`tag'', poly `k`tag'', `inv`tag''/`mass`tag''"
		if `islog`tag'' & `est`tag'' != 4 local lab`tag' "`lab`tag'', log"
	}

	/* total count, to rescale each Dirichlet draw -- always the superset  *
	 * histogram's total, since that is what gets resampled below.        */
	local Ntot = 0
	forvalues i = 1/`nbig' {
		local Ntot = `Ntot' + `Hbig'[`i',1]
	}

	/* ------------------------------------------------------------------ *
	 *  Resample loop (in a private frame so the caller's data is safe)   *
	 * ------------------------------------------------------------------ */
	tempname fr draws
	matrix `draws' = J(`reps',2,.)
	frame create `fr'
	frame `fr' {
		quietly {
			set obs `nbig'
			gen double _y0 = .
			gen double _z  = .
			forvalues i = 1/`nbig' {
				replace _y0 = `Hbig'[`i',1] in `i'
				replace _z  = `Hbig'[`i',2] in `i'
			}
			gen double _y = .
		}

		if "`nodots'" == "" nois _dots 0, ///
			title("Bootstrap replications (refitting both specifications)...") reps(`reps')

		forvalues r = 1/`reps' {
			quietly {
				replace _y = rgamma(_y0, 1)
				summarize _y, meanonly
				replace _y = _y * `Ntot' / r(sum)

				local _ea = .
				capture polbunch _y _z `restrictA', `specA'
				if !_rc {
					capture local _ea = _b[`target']
				}
				if !missing(`_ea') matrix `draws'[`r',1] = `_ea'

				local _eb = .
				capture polbunch _y _z `restrictB', `specB'
				if !_rc {
					capture local _eb = _b[`target']
				}
				if !missing(`_eb') matrix `draws'[`r',2] = `_eb'
			}
			if "`nodots'" == "" nois _dots `r' 0
		}
	}
	frame drop `fr'

	/* ------------------------------------------------------------------ *
	 *  Aggregate                                                          *
	 * ------------------------------------------------------------------ */
	local eA = `_elA'
	local eB = `_elB'
	local d0 = `eA' - `eB'

	tempname fr2
	frame create `fr2'
	local nused = .
	local se    = .
	local corr  = .
	local ppct  = .
	local lb    = .
	local ub    = .
	frame `fr2' {
		quietly {
			set obs `reps'
			svmat double `draws', names(_e)
			gen double _dd = _e1 - _e2
			count if !missing(_dd)
			local nused = r(N)
			if `nused' >= 2 {
				summarize _dd
				local se = r(sd)
				correlate _e1 _e2
				local corr = r(rho)
				local a2 = (100 - `level')/2
				_pctile _dd, p(`a2' `=100-`a2'')
				local lb = r(r1)
				local ub = r(r2)
				count if !missing(_dd) & abs(_dd - `d0') >= abs(`d0')
				local ppct = r(N)/`nused'
			}
		}
	}
	frame drop `fr2'

	local z = .
	local p = .
	if !missing(`se') & `se' > 0 {
		local z = `d0'/`se'
		local p = 2*normal(-abs(`z'))
	}

	/* ------------------------------------------------------------------ *
	 *  Report                                                             *
	 * ------------------------------------------------------------------ */
	di ""
	di as text "Contrast of {res:`tlab'}{txt}: two polbunch estimates on a shared histogram"
	di as text "{hline 74}"
	di as text %-15s abbrev("`nameA'",14) "  " %-36s "`labA'" _col(66) as res %10.4f `eA'
	di as text %-15s abbrev("`nameB'",14) "  " %-36s "`labB'" _col(66) as res %10.4f `eB'
	di as text "{hline 74}"
	di as text %-40s "difference (`=abbrev("`nameA'",14)' - `=abbrev("`nameB'",14)')" _col(66) as res %10.4f `d0'
	if !missing(`se') {
		di as text %-40s "paired bootstrap SE" _col(66) as res %10.4f `se'
		if !missing(`corr') ///
			di as text %-40s "corr across replications" _col(66) as res %10.4f `corr'
		if !missing(`lb') ///
			di as text %-40s "`level'% CI for the difference" ///
			_col(55) as res %8.4f `lb' as text " ," as res %8.4f `ub'
		if !missing(`z') {
			di as text %-40s "H0: difference = 0     z" _col(66) as res %10.4f `z'
			di as text %-40s "p-value (normal / bootstrap)" ///
				_col(66) as res %10.4f `p' as text " /" as res %7.4f `ppct'
		}
		else di as text %-40s "H0: difference = 0     p (bootstrap)" _col(66) as res %10.4f `ppct'
	}
	else {
		di as text "The bootstrap could not form the difference on 2+ replications"
		di as text "(a specification failed to converge on the resamples)."
	}
	di as text "{hline 74}"
	di as text "`nused'/`reps' replications refit both specifications."
	if "`warn'" != "" {
		di as text "Note: the specifications differ in:" as res " `warn'"
		di as text "      the difference reflects those choices as well as the estimator."
	}
	if "`mixnote'" != "" ///
		di as text "Note: `mixnote'."
	di as text "The difference tests whether the two specifications disagree on this sample;"
	di as text "it does not identify which is less biased, and should not be compared with the"
	di as text "difference from a different pair of estimates."

	/* ------------------------------------------------------------------ *
	 *  Return                                                             *
	 * ------------------------------------------------------------------ */
	_estimates unhold `_cur'

	return local target  "`target'"
	return local spec_B  "`labB'"
	return local spec_A  "`labA'"
	return local name_B  "`nameB'"
	return local name_A  "`nameA'"
	if "`mixnote'" != "" return local mixnote "`mixnote'"
	if "`warn'" != ""    return local warning "`warn'"
	return scalar reps       = `reps'
	return scalar reps_used  = `nused'
	if !missing(`ppct') return scalar p_bootstrap = `ppct'
	if !missing(`ub')   return scalar ci_ul = `ub'
	if !missing(`lb')   return scalar ci_ll = `lb'
	if !missing(`corr') return scalar corr  = `corr'
	if !missing(`p')    return scalar p     = `p'
	if !missing(`z')    return scalar z     = `z'
	if !missing(`se')   return scalar se    = `se'
	return scalar value_B = `eB'
	return scalar value_A = `eA'
	return scalar diff    = `d0'
end
