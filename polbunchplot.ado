* polbunchplot version date 20260617
* Author: Martin Eckhoff Andresen
* This program is part of the polbunch package.

capture program drop polbunchplot
program define polbunchplot

    syntax [anything(name=models id="stored estimation name(s)")], ///
        [ graph_opts(string) noci nostar ///
          limit(numlist min=2 max=2) log TRUncate]

    quietly {

        local islog = 0
        capture local islog = e(islog)
        if missing(`islog') local islog = 0

        local nmodels : word count `models'


        /*
        ================================================================
        Zero or one model
        ================================================================
        */

        if `nmodels' == 1 {

            local model : word 1 of `models'

            capture estimates restore `model'
            if _rc {
                local rc = _rc
                noisily display as error ///
                    "Stored estimation result `model' not found."
                exit `rc'
            }
        }


        /*
        ================================================================
        SINGLE MODEL
        ================================================================
        */

        if `nmodels' <= 1 {

            if "`=e(cmdname)'" != "polbunch" {
                noisily display as error ///
                    "Estimates in memory not created by polbunch."
                exit 301
            }

            preserve

            capture confirm matrix e(V)
            if _rc != 0 {
                noisily display as text ///
                    "No variance-covariance matrix found. " ///
                    "Confidence intervals and significance stars not reported."

                local ci   noci
                local star nostar
            }

            clear

            tempname table h0 h1
            tempvar z_est z_orig

            matrix `table' = e(table)

            /*
                IMPORTANT:
                Do NOT use names(col). Some polbunch e(table) matrices
                contain column names that cannot be converted into
                unique Stata variable names.

                names(number) creates table1, table2, ... instead.
            */
            local zcol_first "`e(binname)'"

            local zindex_first = ///
                colnumb(`table', "`zcol_first'")

            if missing(`zindex_first') {
                noisily display as error ///
                    "Could not identify the running-variable column " ///
                    "for first model."
                restore
                exit 498
            }

            /*
                Do not use svmat here. Some Stata versions generate
                temporary names that cannot subsequently be addressed
                reliably. Create the running-variable vector directly
                from the appropriate column of e(table).
            */
            local nrows = rowsof(`table')
            set obs `nrows'

            local zfreq_first = colnumb(`table', "freq")

            generate double `z_est' = .
            generate double freq    = .

            forvalues i = 1/`nrows' {
                replace `z_est' = ///
                    `table'[`i',`zindex_first'] in `i'
                replace freq = ///
                    `table'[`i',`zfreq_first'] in `i'
            }

            if "`e(normalize)'" == "nonormalize" {
                generate double `z_orig' = `z_est'
            }
            else {
                generate double `z_orig' = ///
                    e(cutoff_orig) + ///
                    (`z_est' - e(cutoff_est))*e(xscale)
            }

            local cutoff_plot = e(cutoff_orig)
            local lower_plot  = e(lower_limit)
            local upper_plot  = e(upper_limit)
            local bw_plot     = e(bw_orig)

            if "`e(normalize)'" == "nonormalize" {
                local xarg "x"
            }
            else {
                local cutoff_est_plot = e(cutoff_est)
                local cutoff_org_plot = e(cutoff_orig)
                local xscale_plot     = e(xscale)

                local xarg ///
                    `"(`cutoff_est_plot' + (x - `cutoff_org_plot')/`xscale_plot')"'
            }

            if "`limit'" != "" {

                gettoken min_orig max_orig : limit

                drop if ///
                    `z_orig' < `min_orig' | ///
                    `z_orig' > `max_orig'
            }

            summarize `z_orig', meanonly
            local xmin = r(min)
            local xmax = r(max)

            matrix `h0' = e(b)
            matrix `h0' = `h0'[1,"h0:"]

            matrix `h1' = e(b)
            matrix `h1' = `h1'[1,"h1:"]

            local polynomial = e(polynomial)
            local Kb = `polynomial' + 1

            local h0cons : display %21.17g `h0'[1,`Kb']
            local h1cons : display %21.17g `h1'[1,`Kb']

            local h0cons = strtrim("`h0cons'")
            local h1cons = strtrim("`h1cons'")

            local h0plot `"`h0cons'"'
            local h1plot `"`h1cons'"'

            forvalues i = 1/`polynomial' {

                local h0coef : display %21.17g `h0'[1,`i']
                local h1coef : display %21.17g `h1'[1,`i']

                local h0coef = strtrim("`h0coef'")
                local h1coef = strtrim("`h1coef'")

                local h0plot ///
                    `"`h0plot' + (`h0coef')*(`xarg')^`i'"'

                local h1plot ///
                    `"`h1plot' + (`h1coef')*(`xarg')^`i'"'
            }

            if "`truncate'" != "" {

                summarize freq if ///
                    !inrange(`z_orig', `lower_plot', `upper_plot'), ///
                    meanonly

                replace freq = r(max)*1.2 if freq > r(max)
            }

            summarize freq, meanonly
            local ymax = 1.05*r(max)

            local step = 10^floor(log10(`ymax'/5))

            if `ymax'/`step' > 25 {
                local step = 5*`step'
            }
            else if `ymax'/`step' > 10 {
                local step = 2*`step'
            }

            local ytop = ceil(`ymax'/`step')*`step'

            local yscale ///
                yscale(range(0 `ytop') `log') ///
                ylabel(0(`step')`ytop')

            local zhline

            if `upper_plot' > `cutoff_plot' {
                local zhline ///
                    xline(`upper_plot', ///
                        lcolor(black) ///
                        lpattern(dash))
            }

            local mrline
            local linemax = max(`upper_plot', `xmax')

            capture local MR = ///
                `cutoff_plot' + _b[bunching:marginal_response]

            if _rc == 0 & "`MR'" != "" {

                local mrline ///
                    xline(`MR', ///
                        lcolor(maroon) ///
                        lpattern(longdash))

                local linemax = ///
                    max(`upper_plot', min(`xmax', `MR'))
            }

            if e(estimator) == 4 {

                local h0c = _b[h0:_cons]
                local h1c = _b[h1:_cons]

                local shift = _b[bunching:shift]
                local mr    = _b[bunching:marginal_response]

                local islog = 0
                capture local islog = e(islog)
                if missing(`islog') local islog = 0

                local x0 = `cutoff_plot'
                local x1 = `cutoff_plot' + `mr'

                local y0 = `h0c'

                if `islog' == 1 {
                    local y1 = `h1c'
                }
                else {
                    local y1 = `h1c'/(1 + `shift')
                }

                local linemax = min(`x1', `xmax')

                local hsaez ///
                    `y0' + ///
                    ((`y1' - `y0') / (`x1' - `x0')) * ///
                    (x - `x0')

                twoway ///
                    (bar freq `z_orig', ///
                        barwidth(`bw_plot') ///
                        color(navy%50) ///
                        base(0)) ///
                    (function y=`h0c', ///
                        range(`xmin' `lower_plot') ///
                        lcolor(maroon) ///
                        lpattern(solid)) ///
                    (function y=`h0c', ///
                        range(`lower_plot' `cutoff_plot') ///
                        lcolor(maroon) ///
                        lpattern(shortdash)) ///
                    (function y=`hsaez', ///
                        range(`cutoff_plot' `linemax') ///
                        lcolor(black) ///
                        lpattern(solid)) ///
                    (function y=`h1c', ///
                        range(`cutoff_plot' `upper_plot') ///
                        lcolor(navy) ///
                        lpattern(shortdash)) ///
                    (function y=`h1c', ///
                        range(`upper_plot' `xmax') ///
                        lcolor(navy) ///
                        lpattern(solid)), ///
                    xline(`cutoff_plot', ///
                        lcolor(maroon) ///
                        lpattern(dash)) ///
                    xline(`lower_plot', ///
                        lcolor(black) ///
                        lpattern(dash)) ///
                    xline(`upper_plot', ///
                        lcolor(black) ///
                        lpattern(dash)) ///
                    `mrline' ///
                    graphregion(color(white)) ///
                    plotregion(lcolor(black)) ///
                    ytitle("Frequency") ///
                    xtitle("`zcol'") ///
                    legend( ///
                        label(1 "Frequency") ///
                        label(2 "Estimated h0") ///
                        label(6 "Estimated h1") ///
                        label(4 "Implied counterfactual") ///
                        cols(4) ///
                        order(1 2 6 4) ///
                        pos(6)) ///
                    `yscale' ///
                    `graph_opts'
            }
            else {

                twoway ///
                    (bar freq `z_orig', ///
                        barwidth(`bw_plot') ///
                        color(navy%50) ///
                        base(0)) ///
                    (function y=`h0plot', ///
                        range(`xmin' `lower_plot') ///
                        lcolor(maroon) ///
                        lpattern(solid)) ///
                    (function y=`h0plot', ///
                        range(`lower_plot' `linemax') ///
                        lcolor(maroon) ///
                        lpattern(shortdash)) ///
                    (function y=`h1plot', ///
                        range(`upper_plot' `xmax') ///
                        lcolor(navy) ///
                        lpattern(solid)) ///
                    (function y=`h1plot', ///
                        range(`cutoff_plot' `upper_plot') ///
                        lcolor(navy) ///
                        lpattern(shortdash)), ///
                    xline(`cutoff_plot', ///
                        lcolor(maroon) ///
                        lpattern(dash)) ///
                    xline(`lower_plot', ///
                        lcolor(black) ///
                        lpattern(dash)) ///
                    `mrline' ///
                    `zhline' ///
                    graphregion(color(white)) ///
                    plotregion(lcolor(black)) ///
                    ytitle("Frequency") ///
                    xtitle("`zcol'") ///
                    legend( ///
                        label(1 "Frequency") ///
                        label(2 "Estimated h0") ///
                        label(4 "Estimated h1") ///
                        cols(3) ///
                        order(1 2 4) ///
                        pos(6)) ///
                    `yscale' ///
                    `graph_opts'
            }

            restore
            exit
        }


        /*
        ================================================================
        MULTIPLE MODELS
        ================================================================
        */

        local firstmodel : word 1 of `models'

        capture estimates restore `firstmodel'
        if _rc {
            local rc = _rc
            noisily display as error ///
                "Stored estimation result `firstmodel' not found."
            exit `rc'
        }

        if "`=e(cmdname)'" != "polbunch" {
            noisily display as error ///
                "Stored estimates `firstmodel' were not created by polbunch."
            exit 301
        }

        /*
            There is exactly ONE preserve in the multiple-model branch.
        */
        preserve
        clear


        /*
        ================================================================
        FIRST PASS:
        FIND GLOBAL MINIMUM AND MAXIMUM ACROSS ALL MODELS

        We do this directly from e(table). No svmat, no temporary
        dataset, and therefore no nested preserve/restore.
        ================================================================
        */

        local xmin_global = .
        local xmax_global = .

        foreach model of local models {

            capture estimates restore `model'
            if _rc {
                local rc = _rc
                noisily display as error ///
                    "Stored estimation result `model' not found."
                restore
                exit `rc'
            }

            if "`=e(cmdname)'" != "polbunch" {
                noisily display as error ///
                    "Stored estimates `model' were not created by polbunch."
                restore
                exit 301
            }

            /*
                Copy e(table) to a uniquely named matrix.
            */
            tempname globaltable gmin gmax

            matrix `globaltable' = e(table)

            local zcol_global "`e(binname)'"

            local zindex_global = ///
                colnumb(`globaltable', "`zcol_global'")

            if missing(`zindex_global') {
                noisily display as error ///
                    "Could not identify the running-variable column " ///
                    "for model `model'."
                restore
                exit 498
            }

            /*
                Calculate the min/max in the estimation scale directly
                from the matrix.
            */
            mata: st_numscalar("`gmin'", ///
                min(st_matrix("`globaltable'")[., `zindex_global']))

            mata: st_numscalar("`gmax'", ///
                max(st_matrix("`globaltable'")[., `zindex_global']))

            /*
                Convert to original running-variable scale if needed.
            */
            if "`e(normalize)'" == "nonormalize" {

                local model_min = scalar(`gmin')
                local model_max = scalar(`gmax')
            }
            else {

                local model_min = ///
                    e(cutoff_orig) + ///
                    (scalar(`gmin') - e(cutoff_est))*e(xscale)

                local model_max = ///
                    e(cutoff_orig) + ///
                    (scalar(`gmax') - e(cutoff_est))*e(xscale)
            }

            /*
                Update global range.
            */
            if missing(`xmin_global') {

                local xmin_global = `model_min'
                local xmax_global = `model_max'
            }
            else {

                local xmin_global = ///
                    min(`xmin_global', `model_min')

                local xmax_global = ///
                    max(`xmax_global', `model_max')
            }
        }


        /*
        ================================================================
        SECOND PASS:
        RESTORE FIRST MODEL AND BUILD HISTOGRAM
        ================================================================
        */

        capture estimates restore `firstmodel'
        if _rc {
            local rc = _rc
            noisily display as error ///
                "Stored estimation result `firstmodel' not found."
            restore
            exit `rc'
        }

        tempname table h0 h1
        tempvar z_est z_orig

        matrix `table' = e(table)

        /*
            Do not use svmat here. Create the running-variable vector
            directly from the appropriate column of e(table).
        */
        local zcol_first "`e(binname)'"

        local zindex_first = ///
            colnumb(`table', "`zcol_first'")

        if missing(`zindex_first') {
            noisily display as error ///
                "Could not identify the running-variable column " ///
                "for first model."
            restore
            exit 498
        }

        local zfreq_first = colnumb(`table', "freq")

        local nrows = rowsof(`table')
        set obs `nrows'

        generate double `z_est' = .
        generate double freq    = .

        forvalues i = 1/`nrows' {
            replace `z_est' = ///
                `table'[`i',`zindex_first'] in `i'
            replace freq = ///
                `table'[`i',`zfreq_first'] in `i'
        }

        if "`e(normalize)'" == "nonormalize" {
            generate double `z_orig' = `z_est'
        }
        else {
            generate double `z_orig' = ///
                e(cutoff_orig) + ///
                (`z_est' - e(cutoff_est))*e(xscale)
        }

        local cutoff_first = e(cutoff_orig)
        local lower_first  = e(lower_limit)
        local upper_first  = e(upper_limit)
        local bw_first     = e(bw_orig)


        /*
        ================================================================
        limit()
        ================================================================
        */

        if "`limit'" != "" {

            gettoken min_orig max_orig : limit

            drop if ///
                `z_orig' < `min_orig' | ///
                `z_orig' > `max_orig'

            /*
                The requested limit also restricts the global graph
                range.
            */
            local xmin_global = ///
                max(`xmin_global', `min_orig')

            local xmax_global = ///
                min(`xmax_global', `max_orig')
        }


        /*
        ================================================================
        HISTOGRAM TRUNCATION
        ================================================================
        */

        if "`truncate'" != "" {

            summarize freq if ///
                !inrange(`z_orig', ///
                    `lower_first', ///
                    `upper_first'), ///
                meanonly

            replace freq = r(max)*1.2 if freq > r(max)
        }


        /*
        ================================================================
        Y-AXIS
        ================================================================
        */

        summarize freq, meanonly
        local ymax = 1.05*r(max)

        local step = 10^floor(log10(`ymax'/5))

        if `ymax'/`step' > 25 {
            local step = 5*`step'
        }
        else if `ymax'/`step' > 10 {
            local step = 2*`step'
        }

        local ytop = ceil(`ymax'/`step')*`step'

        local yscale ///
            yscale(range(0 `ytop') `log') ///
            ylabel(0(`step')`ytop')


        /*
        ================================================================
        COLORS
        ================================================================
        */

        local colors ///
            navy maroon forest_green dkorange teal cranberry ///
            purple brown olive sienna ebblue magenta

        local ncolors : word count `colors'


        /*
        ================================================================
        INITIAL PLOT = HISTOGRAM
        ================================================================
        */

        local plots ///
            (bar freq `z_orig', ///
                barwidth(`bw_first') ///
                color(navy%35) ///
                base(0))

        local legend_order 1 "Frequency"

        local plotnum  = 1
        local modelnum = 0


        /*
        ================================================================
        LOOP OVER MODELS
        ================================================================
        */

        foreach model of local models {

            local ++modelnum

            capture estimates restore `model'
            if _rc {
                local rc = _rc
                noisily display as error ///
                    "Stored estimation result `model' not found."
                restore
                exit `rc'
            }

            if "`=e(cmdname)'" != "polbunch" {
                noisily display as error ///
                    "Stored estimates `model' were not created by polbunch."
                restore
                exit 301
            }


            /*
            ============================================================
            GET THIS MODEL'S OWN DATA RANGE

            Again, calculate directly from e(table). There is NO
            preserve/restore here.
            ============================================================
            */

            tempname modeltable mmin mmax

            matrix `modeltable' = e(table)

            local zcol_m "`e(binname)'"

            local zindex_m = ///
                colnumb(`modeltable', "`zcol_m'")

            if missing(`zindex_m') {
                noisily display as error ///
                    "Could not identify the running-variable column " ///
                    "for model `model'."
                restore
                exit 498
            }

            mata: st_numscalar("`mmin'", ///
                min(st_matrix("`modeltable'")[., `zindex_m']))

            mata: st_numscalar("`mmax'", ///
                max(st_matrix("`modeltable'")[., `zindex_m']))

            if "`e(normalize)'" == "nonormalize" {

                local xmin_m = scalar(`mmin')
                local xmax_m = scalar(`mmax')
            }
            else {

                local xmin_m = ///
                    e(cutoff_orig) + ///
                    (scalar(`mmin') - e(cutoff_est))*e(xscale)

                local xmax_m = ///
                    e(cutoff_orig) + ///
                    (scalar(`mmax') - e(cutoff_est))*e(xscale)
            }


            /*
                Apply limit() to this model's range.
            */
            if "`limit'" != "" {

                local xmin_m = ///
                    max(`xmin_m', `min_orig')

                local xmax_m = ///
                    min(`xmax_m', `max_orig')
            }


            /*
            ============================================================
            GET COEFFICIENTS
            ============================================================
            */

            tempname bmodel h0 h1

            matrix `bmodel' = e(b)

            capture matrix `h0' = `bmodel'[1,"h0:"]
            if _rc {
                noisily display as error ///
                    "Stored estimates `model' do not contain h0 coefficients."
                restore
                exit 498
            }

            capture matrix `h1' = `bmodel'[1,"h1:"]
            if _rc {
                noisily display as error ///
                    "Stored estimates `model' do not contain h1 coefficients."
                restore
                exit 498
            }


            /*
            ============================================================
            MODEL PARAMETERS
            ============================================================
            */

            local polynomial_m = e(polynomial)
            local Kb_m         = `polynomial_m' + 1

            local cutoff_m = e(cutoff_orig)
            local lower_m  = e(lower_limit)
            local upper_m  = e(upper_limit)


            /*
            ============================================================
            TRANSFORMATION FOR NORMALIZED ESTIMATION
            ============================================================
            */

            if "`e(normalize)'" == "nonormalize" {

                local xarg_m "x"
            }
            else {

                local cutoff_est_m = e(cutoff_est)
                local cutoff_org_m = e(cutoff_orig)
                local xscale_m     = e(xscale)

                local xarg_m ///
                    `"(`cutoff_est_m' + (x - `cutoff_org_m')/`xscale_m')"'
            }


            /*
            ============================================================
            BUILD H0 AND H1 FUNCTIONS
            ============================================================
            */

            local h0cons_m : display %21.17g ///
                `h0'[1,`Kb_m']

            local h1cons_m : display %21.17g ///
                `h1'[1,`Kb_m']

            local h0cons_m = strtrim("`h0cons_m'")
            local h1cons_m = strtrim("`h1cons_m'")

            local h0plot_m `"`h0cons_m'"'
            local h1plot_m `"`h1cons_m'"'


            forvalues j = 1/`polynomial_m' {

                local h0coef_m : display %21.17g ///
                    `h0'[1,`j']

                local h1coef_m : display %21.17g ///
                    `h1'[1,`j']

                local h0coef_m = strtrim("`h0coef_m'")
                local h1coef_m = strtrim("`h1coef_m'")

                local h0plot_m ///
                    `"`h0plot_m' + (`h0coef_m')*(`xarg_m')^`j'"'

                local h1plot_m ///
                    `"`h1plot_m' + (`h1coef_m')*(`xarg_m')^`j'"'
            }


            /*
            ============================================================
            MODEL COLOR
            ============================================================
            */

            local colorpos = ///
                mod(`modelnum' - 1, `ncolors') + 1

            local modelcolor : word `colorpos' of `colors'


			/*
			================================================================
			MODEL-SPECIFIC PLOTTING RANGES

			First-stage bunching-mass representation:

				h0 solid:
					model minimum -> lower excluded limit

				h0 dotted:
					lower excluded limit -> cutoff

				h1 dotted:
					cutoff -> upper excluded limit

				h1 solid:
					upper excluded limit -> model maximum

			Each segment is restricted to THIS MODEL'S own estimation
			range.
			================================================================
			*/

			local range_h0solid_lo = ///
				`xmin_m'

			local range_h0solid_hi = ///
				min(`lower_m', `xmax_m')


			local range_h0dot_lo = ///
				max(`lower_m', `xmin_m')

			local range_h0dot_hi = ///
				min(`cutoff_m', `xmax_m')


			local range_h1dot_lo = ///
				max(`cutoff_m', `xmin_m')

			local range_h1dot_hi = ///
				min(`upper_m', `xmax_m')


			local range_h1solid_lo = ///
				max(`upper_m', `xmin_m')

			local range_h1solid_hi = ///
				`xmax_m'

            /*
            ------------------------------------------------------------
            h0 below excluded region
            ------------------------------------------------------------
            */

            if `range_h0solid_hi' > `range_h0solid_lo' {

                local ++plotnum
                local legendplot = `plotnum'

                local plots `plots' ///
                    (function y=`h0plot_m', ///
                        range( ///
                            `range_h0solid_lo' ///
                            `range_h0solid_hi' ///
                        ) ///
                        lcolor(`modelcolor') ///
                        lpattern(solid))
            }


            /*
            ------------------------------------------------------------
            h0 THROUGH ENTIRE EXCLUDED REGION
            ------------------------------------------------------------
            */

            if `range_h0dot_hi' > `range_h0dot_lo' {

                local ++plotnum

                if missing(`legendplot') {
                    local legendplot = `plotnum'
                }

                local plots `plots' ///
                    (function y=`h0plot_m', ///
                        range( ///
                            `range_h0dot_lo' ///
                            `range_h0dot_hi' ///
                        ) ///
                        lcolor(`modelcolor') ///
                        lpattern(shortdash))
            }


			/*
------------------------------------------------------------
h0 BELOW excluded region
------------------------------------------------------------
*/

if `range_h0solid_hi' > `range_h0solid_lo' {

    local ++plotnum
    local legendplot = `plotnum'

    local plots `plots' ///
        (function y=`h0plot_m', ///
            range(`range_h0solid_lo' `range_h0solid_hi') ///
            lcolor(`modelcolor') ///
            lpattern(solid))
}


/*
------------------------------------------------------------
h0: LOWER LIMIT -> CUTOFF

This is the h0 reference density used to identify the
bunching mass on the left side of the threshold.
------------------------------------------------------------
*/

if `range_h0dot_hi' > `range_h0dot_lo' {

    local ++plotnum

    if missing(`legendplot') {
        local legendplot = `plotnum'
    }

    local plots `plots' ///
        (function y=`h0plot_m', ///
            range(`range_h0dot_lo' `range_h0dot_hi') ///
            lcolor(`modelcolor') ///
            lpattern(shortdash))
}


	/*
	------------------------------------------------------------
	h1: CUTOFF -> UPPER LIMIT

	This is the h1 reference density used to identify the
	bunching mass on the right side of the threshold.

	IMPORTANT:
	We are NOT using the Saez trapezoid here. This is purely
	the first-stage/reference-density representation.
	------------------------------------------------------------
	*/

	if `range_h1dot_hi' > `range_h1dot_lo' {

		local ++plotnum

		if missing(`legendplot') {
			local legendplot = `plotnum'
		}

		local plots `plots' ///
			(function y=`h1plot_m', ///
				range(`range_h1dot_lo' `range_h1dot_hi') ///
				lcolor(`modelcolor') ///
				lpattern(shortdash))
	}


	/*
	------------------------------------------------------------
	h1 ABOVE excluded region
	------------------------------------------------------------
	*/

	if `range_h1solid_hi' > `range_h1solid_lo' {

		local ++plotnum

		if missing(`legendplot') {
			local legendplot = `plotnum'
		}

		local plots `plots' ///
			(function y=`h1plot_m', ///
				range(`range_h1solid_lo' `range_h1solid_hi') ///
				lcolor(`modelcolor') ///
				lpattern(solid))
	}

            /*
            ------------------------------------------------------------
            h1 STARTS AT UPPER EXCLUDED LIMIT
            ------------------------------------------------------------
            */

            if `range_h1solid_hi' > `range_h1solid_lo' {

                local ++plotnum

                if missing(`legendplot') {
                    local legendplot = `plotnum'
                }

                local plots `plots' ///
                    (function y=`h1plot_m', ///
                        range( ///
                            `range_h1solid_lo' ///
                            `range_h1solid_hi' ///
                        ) ///
                        lcolor(`modelcolor') ///
                        lpattern(solid))
            }


            if missing(`legendplot') {
                noisily display as error ///
                    "No part of model `model' lies inside the plotting range."
                restore
                exit 498
            }


            /*
            ============================================================
            LEGEND
            ============================================================
            */

            local model_label = ///
                upper(substr("`model'",1,1)) + ///
                substr("`model'",2,.)

            local legend_order ///
                `legend_order' ///
                `legendplot' ///
                "`model_label'"
        }


        /*
        ================================================================
        REFERENCE LINES
        ================================================================
        */

        local lowerline

        if `lower_first' < `cutoff_first' {

            local lowerline ///
                xline(`lower_first', ///
                    lcolor(black) ///
                    lpattern(dash))
        }

        local upperline

        if `upper_first' > `cutoff_first' {

            local upperline ///
                xline(`upper_first', ///
                    lcolor(black) ///
                    lpattern(dash))
        }


        local legend_cols = min(`nmodels' + 1, 5)


        /*
        ================================================================
        FINAL GRAPH

        Histogram/x-axis:
            GLOBAL min -> GLOBAL max across models.

        Polynomial curves:
            each model's OWN estimation range.
        ================================================================
        */

        twoway ///
            `plots', ///
            xscale(range(`xmin_global' `xmax_global')) ///
            xline(`cutoff_first', ///
                lcolor(maroon) ///
                lpattern(dash)) ///
            `lowerline' ///
            `upperline' ///
            graphregion(color(white)) ///
            plotregion(lcolor(black)) ///
            ytitle("Frequency") ///
            xtitle("`zcol_first'") ///
            legend( ///
                order(`legend_order') ///
                cols(`legend_cols') ///
                pos(6)) ///
            `yscale' ///
            `graph_opts'

        /*
            THE ONLY restore in the multiple-model branch.
        */
        restore
    }
end