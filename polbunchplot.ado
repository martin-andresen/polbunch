*! polbunchplot version date 20241210
* Author: Martin Eckhoff Andresen
* This program is part of the polbunch package.
cap prog drop polbunchplot
	program polbunchplot
	
	syntax [name], [graph_opts(string) noci nostar limit(numlist min=2 max=2) log TRUncate]
	quietly {
		if "`name'"!="" est restore `name'
		if "`=e(cmdname)'"!="polbunch" {
			noi di in red "Estimates in memory not created by polbunch."
			exit
			}
		
		preserve
		
		cap confirm matrix e(V)
		if _rc!=0 {
			noi di as text "No variance-covariance matrix found. Confidence intervals and significance stars not reported."
			loc ci noci
			loc star nostar
		}
		
		clear
		
		tempvar h0 h1 freq
		mat `freq'= e(table)
		svmat `freq', names(col)
		
		if "`limit'"!="" {
			gettoken min max: limit
			drop if `e(binname)'<`min'|`e(binname)'>`max'
		}
		su `e(binname)'
		loc xmin=r(min)
		loc xmax=r(max)
		
		loc plus=0
		mat `h0'=e(b)
		mat `h0'=`h0'[1,"h0:"]
		mat `h1'=e(b)
		mat `h1'=`h1'[1,"h1:"]
		
		tempname h0plot h1plot
		loc h0plot `h0'[1,`=colsof(`h0')']
		loc h1plot `h1'[1,`=colsof(`h1')']
		forvalues i=1/`=e(polynomial)' {
			loc h0plot `h0plot' +`h0'[1,`i']*x^`i' 
			loc h1plot `h1plot' +`h1'[1,`i']*x^`i' 
		}
		
		if "`truncate'"!="" {
			su freq if !inrange(`=e(binname)',`=e(lower_limit)',`=e(upper_limit)')
			replace freq=r(max) if freq>r(max)
		}
		if "`log'"=="log" loc yscale yscale(log)
		if e(upper_limit)>e(cutoff) loc zhline xline(`=e(upper_limit)', lcolor(black) lpattern(dash))
		twoway 	(bar freq `=e(binname)', barwidth(`=e(bw)') color(navy%50)) ///
				(function h0=`h0plot', range(`xmin' `=e(lower_limit)') lcolor(maroon) lpattern(solid)) ///
				(function h0=`h0plot', range(`=e(lower_limit)' `xmax') lcolor(maroon) lpattern(shortdash)) ///
				(function h1=`h1plot', range(`=e(upper_limit)' `xmax') lcolor(navy) lpattern(solid)), ///
				xline(`=e(cutoff)', lcolor(maroon) lpattern(dash)) xline(`=e(lower_limit)', lcolor(black) lpattern(dash)) `zhline' ///
				graphregion(color(white)) plotregion(lcolor(black)) ytitle("Frequency") xtitle("`=e(binname)'") ///
				legend(label(1 "Frequency") label(2 "Estimated h0") label(4 "Estimated h1") cols(3) order(1 2 4) pos(6)) `yscale'
		
		restore
	}
	
	end
