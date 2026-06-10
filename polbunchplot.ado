*! polbunchplot version date 20241210
* Author: Martin Eckhoff Andresen
* This program is part of the polbunch package.

cap prog drop polbunchplot
program polbunchplot

    syntax [name], [graph_opts(string) noci nostar limit(numlist min=2 max=2) log TRUncate]

    quietly {
        if "`name'" != "" est restore `name'

        if "`=e(cmdname)'" != "polbunch" {
            noi di in red "Estimates in memory not created by polbunch."
            exit
        }

        preserve

        cap confirm matrix e(V)
        if _rc != 0 {
            noi di as text "No variance-covariance matrix found. Confidence intervals and significance stars not reported."
            local ci noci
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
                z in estimation scale

            Rename the stored z column to z_est and create z_orig for plotting.
        */
        local zcol "`e(binname)'"
        rename `zcol' `z_est'

        if "`e(normalize)'" == "nonormalize" {
            gen double `z_orig' = `z_est'

            local cutoff_plot = e(cutoff)
            local lower_plot  = e(lower_limit)
            local upper_plot  = e(upper_limit)
            local bw_plot     = e(bw)

            local xarg "x"
        }
        else {
            gen double `z_orig' = e(cutoff) + e(bw)*`z_est'

            local cutoff_plot = e(cutoff)
            local lower_plot  = e(lower_limit)
            local upper_plot  = e(upper_limit)
            local bw_plot     = e(bw)

            local xarg "((x - e(cutoff))/e(bw))"
        }

        /*
            Interpret limit() in original scale.
        */
        if "`limit'" != "" {
            gettoken min_orig max_orig : limit
            drop if `z_orig' < `min_orig' | `z_orig' > `max_orig'
        }

        summarize `z_orig', meanonly
        local xmin = r(min)
        local xmax = r(max)

        matrix `h0' = e(b)
        matrix `h0' = `h0'[1,"h0:"]

        matrix `h1' = e(b)
        matrix `h1' = `h1'[1,"h1:"]

        /*
            Coefficients are in estimation scale.
            Functions are plotted over original-scale x, so evaluate them at xarg.
        */
        local h0plot `h0'[1,`=colsof(`h0')']
        local h1plot `h1'[1,`=colsof(`h1')']

        forvalues i = 1/`=e(polynomial)' {
            local h0plot `h0plot' + `h0'[1,`i']*(`xarg')^`i'
            local h1plot `h1plot' + `h1'[1,`i']*(`xarg')^`i'
        }

        if "`truncate'" != "" {
            summarize freq if !inrange(`z_orig', `lower_plot', `upper_plot'), meanonly
            replace freq = r(max) if freq > r(max)
        }

		su freq, meanonly
		local ymax = 1.05*r(max)

		local step = 10^floor(log10(`ymax'/5))
		if `ymax'/`step' > 25 local step = 5*`step'
		else if `ymax'/`step' > 10 local step = 2*`step'

		local ytop = ceil(`ymax'/`step')*`step'

		local yscale yscale(range(0 `ytop') `log') ylabel(0(`step')`ytop')

        if `upper_plot' > `cutoff_plot' {
            local zhline xline(`upper_plot', lcolor(black) lpattern(dash))
        }

		cap loc MR=`cutoff_plot'+_b[bunching:marginal_response]
		if !mi(`MR') loc mrline xline(`MR', lcolor(maroon) lpattern(longdash))
        if e(estimator)==4 {
			local h0c = _b[h0:_cons]
			local h1c = _b[h1:_cons]

			local htrap `h0c' + ///
				((`h1c' - `h0c') / (`upper_plot' - `lower_plot')) * ///
				(x - `lower_plot')

			twoway ///
				(bar freq `z_orig', barwidth(`bw_plot') color(navy%50) base(0)) ///
				(function y=`h0c', range(`xmin' `lower_plot') lcolor(maroon) lpattern(solid)) ///
				(function y=`htrap', range(`lower_plot' `upper_plot') lcolor(maroon) lpattern(shortdash)) ///
				(function y=`h1c', range(`upper_plot' `xmax') lcolor(navy) lpattern(solid)), ///
				xline(`cutoff_plot', lcolor(maroon) lpattern(dash)) ///
				xline(`lower_plot',  lcolor(black)  lpattern(dash)) ///
				xline(`upper_plot',  lcolor(black)  lpattern(dash)) ///
				graphregion(color(white)) ///
				plotregion(lcolor(black)) ///
				ytitle("Frequency") ///
				xtitle("`zcol'") ///
				`mrline' ///
				legend(label(1 "Frequency") label(2 "Estimated h0") label(4 "Estimated h1") ///
                   cols(3) order(1 2 4) pos(6)) ///
				`yscale' ///
				`graph_opts'
		}
		else {
			twoway ///
            (bar freq `z_orig', barwidth(`bw_plot') color(navy%50) base(0)) ///
            (function h0=`h0plot', range(`xmin' `lower_plot') lcolor(maroon) lpattern(solid)) ///
            (function h0=`h0plot', range(`lower_plot' `MR') lcolor(maroon) lpattern(shortdash)) ///
            (function h1=`h1plot', range(`upper_plot' `xmax') lcolor(navy) lpattern(solid)) ///
			(function h1=`h1plot', range(`cutoff_plot' `upper_plot') lcolor(navy) lpattern(shortdash)), ///
            xline(`cutoff_plot', lcolor(maroon) lpattern(dash)) ///
            xline(`lower_plot', lcolor(black) lpattern(dash)) ///
			`mrline' ///
            `zhline' ///
            graphregion(color(white)) ///
            plotregion(lcolor(black)) ///
            ytitle("Frequency") ///
            xtitle("`zcol'") ///
            legend(label(1 "Frequency") label(2 "Estimated h0") label(4 "Estimated h1") ///
                   cols(3) order(1 2 4) pos(6)) ///
            `yscale' ///
            `graph_opts'
		}
        restore
    }

end