*! polbunchplot version date 20260617
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
            No supplied name:
                Use estimates currently in memory.

            One supplied name:
                Restore that stored estimation result and use the
                single-model plotting routine.

            Multiple supplied names:
                Use the first model for the histogram and overlay
                one composite h0/h1 curve for every stored model.
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
        Zero or one model: single-model plot
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
            svmat `table', names(col)

            /*
                e(table) contains:
                    freq
                    running variable in estimation scale
            */
            local zcol "`e(binname)'"
            rename `zcol' `z_est'

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

            /*
                Map original-scale graph coordinate x into the
                estimation scale used for the polynomial coefficients.
            */
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

            /*
                limit() is interpreted in the original running-variable
                scale.
            */
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

            /*
                Convert coefficients to numeric text and build graph
                function expressions without leaving matrix references
                or literal quotation marks in the expressions.
            */
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
					xline(`cutoff_plot', lcolor(maroon) lpattern(dash)) ///
					xline(`lower_plot', lcolor(black) lpattern(dash)) ///
					xline(`upper_plot', lcolor(black) lpattern(dash)) ///
					`mrline' ///
					graphregion(color(white)) ///
					plotregion(lcolor(black)) ///
					ytitle("Frequency") ///
					xtitle("`zcol'") ///
					legend(label(1 "Frequency") label(2 "Estimated h0") ///
						   label(6 "Estimated h1") label(4 "Implied counterfactual") cols(4) order(1 2 6 4) pos(6)) ///
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
                        pos(6) ///
                    ) ///
                    `yscale' ///
                    `graph_opts'
            }

            restore
            exit
        }


        /*
        ================================================================
        Multiple models
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

        preserve
        clear

        tempname table bmodel h0 h1
        tempvar z_est z_orig

        /*
            The first model supplies the histogram and plotting sample.
        */
        matrix `table' = e(table)
        svmat `table', names(col)

        local zcol_first "`e(binname)'"
        rename `zcol_first' `z_est'

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
            limit() is interpreted in the original scale.
        */
        if "`limit'" != "" {

            gettoken min_orig max_orig : limit

            drop if ///
                `z_orig' < `min_orig' | ///
                `z_orig' > `max_orig'
        }

        summarize `z_orig', meanonly
        local xmin = r(min)
        local xmax = r(max)

        if "`truncate'" != "" {

            summarize freq if ///
                !inrange(`z_orig', `lower_first', `upper_first'), ///
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

        /*
            Colors are recycled if more models are supplied.
        */
        local colors ///
            navy maroon forest_green dkorange teal cranberry ///
            purple brown olive sienna ebblue magenta

        local ncolors : word count `colors'

        /*
            Plot 1 is the histogram.
        */
        local plots ///
            (bar freq `z_orig', ///
                barwidth(`bw_first') ///
                color(navy%35) ///
                base(0))

        local legend_order 1 "Frequency"
        local plotnum  = 1
        local modelnum = 0

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

            local polynomial_m = e(polynomial)
            local Kb_m         = `polynomial_m' + 1

            local cutoff_m = e(cutoff_orig)
            local lower_m  = e(lower_limit)
            local upper_m  = e(upper_limit)

            /*
                Map original-scale graph coordinate x into this model's
                estimation scale.
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
                Convert coefficients to trimmed numeric literals and
                construct the model's h0 and h1 expressions.
            */
            local h0cons_m : display %21.17g `h0'[1,`Kb_m']
            local h1cons_m : display %21.17g `h1'[1,`Kb_m']

            local h0cons_m = strtrim("`h0cons_m'")
            local h1cons_m = strtrim("`h1cons_m'")

            local h0plot_m `"`h0cons_m'"'
            local h1plot_m `"`h1cons_m'"'

            forvalues j = 1/`polynomial_m' {

                local h0coef_m : display %21.17g `h0'[1,`j']
                local h1coef_m : display %21.17g `h1'[1,`j']

                local h0coef_m = strtrim("`h0coef_m'")
                local h1coef_m = strtrim("`h1coef_m'")

                local h0plot_m ///
                    `"`h0plot_m' + (`h0coef_m')*(`xarg_m')^`j'"'

                local h1plot_m ///
                    `"`h1plot_m' + (`h1coef_m')*(`xarg_m')^`j'"'
            }

            local colorpos = ///
                mod(`modelnum' - 1, `ncolors') + 1

            local modelcolor : word `colorpos' of `colors'

            /*
                One composite model curve:

                    h0 solid below the excluded region;
                    h0 dotted from lower excluded limit to cutoff;
                    h1 dotted from cutoff to upper excluded limit;
                    h1 solid above the excluded region.

                Only the first available segment gets a legend entry.
            */
            local range_h0solid_lo = `xmin'
            local range_h0solid_hi = min(`lower_m', `xmax')

            local range_h0dot_lo = max(`lower_m', `xmin')
            local range_h0dot_hi = min(`cutoff_m', `xmax')

            local range_h1dot_lo = max(`cutoff_m', `xmin')
            local range_h1dot_hi = min(`upper_m', `xmax')

            local range_h1solid_lo = max(`upper_m', `xmin')
            local range_h1solid_hi = `xmax'

            local legendplot = .

            if `range_h0solid_hi' > `range_h0solid_lo' {

                local ++plotnum
                local legendplot = `plotnum'

                local plots `plots' ///
                    (function y=`h0plot_m', ///
                        range(`range_h0solid_lo' `range_h0solid_hi') ///
                        lcolor(`modelcolor') ///
                        lpattern(solid))
            }

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

            if missing(`legendplot') {
                noisily display as error ///
                    "No part of model `model' lies inside the plotting range."
                restore
                exit 498
            }

            local model_label = upper(substr("`model'",1,1)) + substr("`model'",2,.)
            local legend_order ///
                `legend_order' `legendplot' "`model_label'"
        }

        /*
            Reference lines use the first model, which also supplies
            the histogram.
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

        twoway ///
            `plots', ///
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
                pos(6) ///
            ) ///
            `yscale' ///
            `graph_opts'

        restore
    }
end