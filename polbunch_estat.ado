*! polbunch_estat version date 20260910
* Author: Martin Eckhoff Andresen
* estat post-estimation for polbunch: reference-region goodness of fit.
* Registered via  ereturn local estat_cmd "polbunch_estat"  in polbunch.ado.
*
*   estat gof                       full block for the active polbunch estimate
*   estat gof m1 m2 m3 [, common()] comparison table across stored estimates,
*                                   with a common phi-hat for the QAICc column

cap program drop polbunch_estat
program define polbunch_estat, rclass
	version 16.0

	gettoken sub 0 : 0, parse(" ,")
	local sub = strlower(strtrim("`sub'"))
	if !inlist("`sub'", "gof", "") {
		di as error "estat subcommand `sub' is not available after polbunch (available: gof)"
		exit 198
	}

	syntax [anything(name=models)] [, COMMON(string) ]

	if `"`models'"' != "" {
		_polbunch_estat_gofcmp `models', common(`"`common'"')
		return add
		exit
	}

	*---------------------------------------------------------------
	* single active estimate
	*---------------------------------------------------------------
	if "`e(cmd)'" != "polbunch" {
		di as error "last estimates are not polbunch results"
		exit 301
	}
	if missing(e(r2_dev)) {
		di as error "no goodness-of-fit statistics stored for these results"
		di as error "(vce(bootstrap, multinomial|bayesian) and vce(none) do not compute them)"
		exit 301
	}

	local df = e(gof_df)
	local phi = .
	if `df' > 0 & `df' < . local phi = e(pearson_x2) / `df'

	di as txt _n "Reference-region goodness of fit  (Poisson / quasi-Poisson family)"
	di as txt "{hline 64}"
	di as txt %-36s "Bins fit  /  parameters  /  df" ///
		as res %6.0f e(gof_nbins) "  /" %6.0f e(gof_np) "  /" %6.0f `df'
	di as txt %-36s "Deviance" as res %13.3f e(deviance) ///
		as txt "   p = " as res %5.3f e(deviance_p)
	di as txt %-36s "Pearson X2" as res %13.3f e(pearson_x2) ///
		as txt "   p = " as res %5.3f e(pearson_x2_p)
	di as txt %-36s "Pearson dispersion  (phi-hat)" as res %13.3f `phi'
	di as txt %-36s "Deviance R2  (Cameron-Windmeijer)" as res %13.3f e(r2_dev)
	di as txt %-36s "RMS reference-bin residual" as res %13.3f e(gof_rmse)
	if !missing(e(dispersion_below)) | !missing(e(dispersion_above)) {
		di as txt "{hline 64}"
		di as txt %-36s "  phi-hat below cutoff  (A3 only)" as res %13.3f e(dispersion_below) ///
			as txt "  df " as res %5.0f e(gof_df_below)
		di as txt %-36s "  phi-hat above cutoff  (+ A1/A2)" as res %13.3f e(dispersion_above) ///
			as txt "  df " as res %5.0f e(gof_df_above)
	}
	di as txt "{hline 64}"
	di as txt %-36s "Log-likelihood (Poisson)" as res %13.2f e(gof_ll)
	di as txt %-36s "AIC  /  BIC" as res %13.1f e(aic) as res "  /" %13.1f e(bic)
	di as txt %-36s "QAIC  /  QAICc   (own phi-hat)" as res %13.1f e(qaic) as res "  /" %13.1f e(qaicc)
	di as txt "{hline 64}"
	di as txt "Note: the deviance / Pearson p-values reject at near-census N."
	di as txt "      Compare R2, QAICc and residual structure across models."
	di as txt "      For a formal multi-model QAIC ranking use one common"
	di as txt "      phi-hat --  estat gof <models>, common(<most general model>)"

	return scalar ll               = e(gof_ll)
	return scalar df               = `df'
	return scalar dispersion       = `phi'
	return scalar dispersion_below = e(dispersion_below)
	return scalar dispersion_above = e(dispersion_above)
	return scalar bic        = e(bic)
	return scalar aic        = e(aic)
	return scalar qaicc      = e(qaicc)
	return scalar qaic       = e(qaic)
	return scalar r2_dev     = e(r2_dev)
	return scalar pearson_x2 = e(pearson_x2)
	return scalar deviance   = e(deviance)
end


*=================================================================
* Comparison table across stored polbunch estimates.
*   common(name) : use that model's phi-hat for every QAICc
*   common(#)    : use a literal phi-hat
*   (omitted)    : use the listed model with the most parameters
*                  (Burnham-Anderson: phi-hat from the global model)
* All models must have been fit to the same reference bins.
*=================================================================
cap program drop _polbunch_estat_gofcmp
program define _polbunch_estat_gofcmp, rclass
	version 16.0
	syntax anything(name=models) [, COMMON(string) ]

	tempname hold
	capture _estimates hold `hold', restore nullok

	*--- resolve the common phi-hat -----------------------------------
	local phi_c = .
	local phi_src ""
	if `"`common'"' != "" {
		capture confirm number `common'
		if !_rc {
			local phi_c  = `common'
			local phi_src "supplied"
		}
		else {
			capture estimates restore `common'
			if _rc {
				di as error "common(`common') is neither a number nor a stored estimate"
				capture _estimates unhold `hold'
				exit 198
			}
			local phi_c  = e(dispersion)
			local phi_src "`common'"
		}
	}
	if missing(`phi_c') {
		local bestnp = -1
		foreach m of local models {
			capture estimates restore `m'
			if _rc continue
			if "`e(cmd)'" == "polbunch" & !missing(e(gof_np)) & e(gof_np) > `bestnp' {
				local bestnp  = e(gof_np)
				local phi_c   = e(dispersion)
				local phi_src "`m'"
			}
		}
	}
	if missing(`phi_c') local phi_c = 1
	if `phi_c' < 1      local phi_c = 1

	*--- collect the rows -------------------------------------------
	local r  = 0
	local n0 = .
	foreach m of local models {
		capture estimates restore `m'
		if _rc {
			di as txt "  (skipping `m': not a stored estimate)"
			continue
		}
		if "`e(cmd)'" != "polbunch" | missing(e(gof_ll)) {
			di as txt "  (skipping `m': no polbunch goodness of fit stored)"
			continue
		}
		if `n0' == . local n0 = e(gof_nbins)
		else if e(gof_nbins) != `n0' {
			di as error "  `m': `=e(gof_nbins)' reference bins, expected `n0'"
			di as error "  models fit to different windows are not comparable"
			capture _estimates unhold `hold'
			exit 198
		}
		local ++r
		local rn_`r'    "`m'"
		local est_`r'   = e(estimator)
		local K_`r'     = e(polynomial)
		local np_`r'    = e(gof_np)
		local df_`r'    = e(gof_df)
		local dev_`r'   = e(deviance)
		local r2_`r'    = e(r2_dev)
		local q_`r'     = -2*e(gof_ll)/`phi_c' + 2*e(gof_np)
		local penc      = cond(`n0'-e(gof_np)-1 > 0, 2*e(gof_np)*(e(gof_np)+1)/(`n0'-e(gof_np)-1), 0)
		local qc_`r'    = `q_`r'' + `penc'
	}
	if `r' == 0 {
		di as error "no comparable polbunch estimates in the list"
		capture _estimates unhold `hold'
		exit 498
	}

	local qmin = .
	forvalues i = 1/`r' {
		if `qc_`i'' < `qmin' local qmin = `qc_`i''
	}
	local wsum = 0
	forvalues i = 1/`r' {
		local d_`i' = `qc_`i'' - `qmin'
		local w_`i' = exp(-`d_`i''/2)
		local wsum  = `wsum' + `w_`i''
	}

	*--- print -----------------------------------------------------
	di as txt _n "Reference-region goodness of fit across models"
	di as txt "common phi-hat = " as res %5.2f `phi_c' as txt "  (`phi_src')" ///
		"   reference bins = " as res `n0'
	di as txt "{hline 76}"
	di as txt %-12s "Model" %4s "Est" %4s "K" %7s "df" %13s "Deviance" ///
		%9s "R2_dev" %11s "QAICc" %9s "dQAICc" %9s "weight"
	di as txt "{hline 76}"
	forvalues i = 1/`r' {
		di as res %-12s "`rn_`i''" %4.0f `est_`i'' %4.0f `K_`i'' %7.0f `df_`i'' ///
			%13.2f `dev_`i'' %9.3f `r2_`i'' %11.1f `qc_`i'' %9.1f `d_`i'' ///
			%9.3f `w_`i''/`wsum'
	}
	di as txt "{hline 76}"
	di as txt "dQAICc from the best model in the list; weight = Akaike weight."
	di as txt "QAICc uses the common phi-hat above, not each model's own."

	*--- return --------------------------------------------------
	tempname M
	matrix `M' = J(`r', 8, .)
	local rown ""
	forvalues i = 1/`r' {
		matrix `M'[`i',1] = `est_`i''
		matrix `M'[`i',2] = `K_`i''
		matrix `M'[`i',3] = `df_`i''
		matrix `M'[`i',4] = `dev_`i''
		matrix `M'[`i',5] = `r2_`i''
		matrix `M'[`i',6] = `qc_`i''
		matrix `M'[`i',7] = `d_`i''
		matrix `M'[`i',8] = `w_`i''/`wsum'
		local rown "`rown' `rn_`i''"
	}
	matrix colnames `M' = estimator K df deviance r2_dev qaicc dqaicc weight
	matrix rownames `M' = `rown'
	return matrix table = `M'
	return scalar phi_common = `phi_c'
	return local  common_src "`phi_src'"

	capture _estimates unhold `hold'
end
