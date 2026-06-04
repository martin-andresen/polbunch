// bunchcalc: Transforms estimates to h0,h1,B,excess_mass,shift,marginal_response,(elasticity)
program bunchcalc, rclass
	syntax, estimator(integer) cutoff(real) polynomial(integer) bw(real) h(integer) l(integer) [b0(string) 	t0(numlist min=1 max=1 <1) t1(numlist min=1 max=1 <1) constant nopositiveshift boot log nlcom(string) nonormalize]
	
	if "`positiveshift'"!="nopositiveshift" loc shift exp(_b[/lnshift])
	else loc shift _b[/shift]
	
	// cutoff() and bw() are always passed in original units.
	// If normalized, the fitted polynomial is in x = (z-cutoff)/bw units.
	// Then the estimation-scale cutoff/bw are 0/1, and xscale converts
	// marginal response back to original z units.
	local cutoff_orig = `cutoff'
	local bw_orig     = `bw'

	if "`normalize'" == "nonormalize" {
		local cutoff_est = `cutoff_orig'
		local bw_est     = `bw_orig'
		local xscale     = 1
	}
	else {
		local cutoff_est = 0
		local bw_est     = 1
		local xscale     = `bw_orig'
	}

	//b0 (constant in h0)
	if "`b0'"!=""&`estimator'>0 {
		if "`positiveshift'"!="nopositiveshift" loc b0=subinstr("`b0'","{lnshift}","_b[/lnshift]",.)
		else loc b0=subinstr("`b0'","{shift}","_b[/shift]",.)
		forvalues k=1/`=`polynomial'' {
			loc b0=subinstr("`b0'","{b`k'}","_b[/b`k']",.)
		}
		forvalues b=1/`=`h'+`l'' {
			loc b0=subinstr("`b0'","{bunch`b'}","_b[/bunch`b']",.)
		}
	}
	
	//number bunchers
	loc B=0
	forvalues b=1/`=`h'+`l'' {
		loc B `B'+_b[/bunch`b']
		}
	
	//compute nlcom string
	if "`nlcom'"=="" {
		//bk (parameters in h0)
		forvalues k=1/`=`polynomial'' {
			loc nlcom `nlcom' (_b[/b`k'])
		}
		loc nlcom `nlcom' (`b0')
		
		//gk (parameters in h1)
		forvalues k=1/`polynomial' {
			if `estimator'==0 loc nlcom `nlcom' (_b[/g`k'])
			else if `estimator'==1 loc nlcom `nlcom' (_b[/b`k'])
			else if `estimator'==2 loc nlcom `nlcom' (_b[/b`k']/(1+`shift'))
			else if "`log'"=="" loc nlcom `nlcom' (_b[/b`k']*(1+`shift')^(`=`k'+1'))
			else {
					loc str _b[/b`k']
					if (`polynomial'>`k') {
						forvalues n=`=`k'+1'/`=`polynomial'' {
							loc str `str' +_b[/b`n']*comb(`n',`k')*ln(1+`shift')^(`n'-`k')
						}
					}
				loc nlcom `nlcom' (`str')
				}
			}
		
		//g0 (constant in h1)
		if `estimator'==0 loc nlcom `nlcom' (_b[/g0])
		else if `estimator'==1 loc nlcom `nlcom' (`b0')
		else if `estimator'==2 loc nlcom `nlcom' ((`b0')/(1+`shift'))
		else if "`log'"=="" loc nlcom `nlcom' ((`b0')*(1+`shift')) 
		else {
			loc str `b0'
			forvalues n=1/`=`polynomial'' {
					loc str `str' +_b[/b`n']*ln(1+`shift')^`n'
				}
				loc nlcom `nlcom' (`str')
		}
		
		//delta from chetty estimator
		if `estimator'==2 {
			loc nlcom `nlcom' (`shift')
		}

		//number bunchers
		loc nlcom `nlcom' (`B')
		loc m `b0'
		forvalues k=1/`polynomial' {
			loc m `m' + _b[/b`k']*`cutoff_est'^(`k')
		}
		loc nlcom `nlcom' ((`B')/(`m')) //excess mass
		
		if "`constant'"!=""	loc deltaz (`bw_orig'*(`B')/(`m'))
		
		if `estimator'==3&"`constant'"=="" {
			loc deltaz `shift'*`cutoff_orig'
			loc nlcom `nlcom' (`shift') (`deltaz')
			if "`t0'"!=""&"`t1'"!=""  loc nlcom `nlcom' (ln(1+`shift')/(ln(1-(`t0'))-ln(1-(`t1'))))
			}
		else if "`constant'"!="" {
			if "`log'"=="" loc nlcom `nlcom' (`deltaz'/`cutoff_orig') (`deltaz')
			else loc nlcom `nlcom' (exp(`deltaz'+`cutoff')/exp(`cutoff_orig')) (`deltaz')
			if "`t0'"!=""&"`t1'"!="" { //elasticity
				if "`log'"=="" loc nlcom `nlcom' (ln(1+ `deltaz'/`cutoff_orig')/(ln(1-(`t0'))-ln(1-(`t1')))) 
				else loc nlcom `nlcom' ((ln(exp(`cutoff'+`deltaz'))-`cutoff_orig')/(ln(1-(`t0'))-ln(1-(`t1'))))
			}
		}
	}
			
			
	/// TRANSFORM ESTIMATES
	tempname b
	if "`boot'"=="" {
		tempname V
		nlcom `nlcom'
		mat `b'=r(b)
		mat `V'=r(V)
		if `estimator'<3&"`constant'"=="" {
			//calculate shift/MR/el without variance and add to b,V
			tempname h0
			mat `h0'=`b0'
			forvalues k=1/`polynomial' {
				mat `h0'=`h0',_b[/b`k']
			}
			mata: st_numscalar("eresp",eresp(`=`B'',`cutoff_est',st_matrix("`h0'"),`bw_est',`xscale'))
			if eresp==. loc exit=1
			
			if "`log'"=="" mat `b'=`b',`=eresp/`cutoff_orig'',eresp
			else mat `b'=`b',`=exp(eresp+`cutoff')/exp(`cutoff')',eresp
			if "`t0'"!=""&"`t1'"!="" { //elasticity
				if "`log'"=="" mat `b'=`b',`=ln(1+ eresp/`cutoff_orig')/(ln(1-(`t0'))-ln(1-(`t1')))'
				else mat `b'=`b',`=(ln(exp(`cutoff'+eresp))-`cutoff_orig')/(ln(1-(`t0'))-ln(1-(`t1')))'
			}

			if "`t0'"!=""&"`t1'"!="" loc extra=3
			else loc extra=2
			mat `V'=[`V', J(rowsof(`V'),`extra',0) \ J(`extra',colsof(`V'),0) , J(`extra',`extra',0)]
		}
	}
	else {
		 loc nlcom2 `nlcom'
		 while "`nlcom2'"!="" {
			gettoken use nlcom2: nlcom2, match(parns)
			mat `b'=nullmat(`b'),`=`use''
		}
		if `estimator'<3&"`constant'"=="" {
			tempname h0
			mat `h0'=`b0'
			forvalues k=1/`polynomial' {
				mat `h0'=`h0',_b[/b`k']
			}
			mata: st_numscalar("eresp",eresp(`=`B'',`cutoff_est',st_matrix("`h0'"),`bw_est',`xscale'))
			if eresp==. loc exit=1
			if "`log'"=="" mat `b'=`b',`=eresp/`cutoff_orig'',eresp
			else mat `b'=`b',`=exp(eresp+`cutoff')/exp(`cutoff')',eresp
			if "`t0'"!=""&"`t1'"!="" { //elasticity
				if "`log'"=="" mat `b'=`b',`=ln(1+ eresp/`cutoff_orig')/(ln(1-(`t0'))-ln(1-(`t1')))'
				else mat `b'=`b',`=(ln(exp(`cutoff'+eresp))-`cutoff_orig')/(ln(1-(`t0'))-ln(1-(`t1')))'
			}

		}
	}
	
	//NAMES!!
	return local nlcom "`nlcom'"
	return matrix b=`b'
	if "`exit'"=="" loc exit=0
	return scalar exit=`exit'
	if "`boot'"=="" {
		return matrix V=`V'
	}
	
	end
	

mata:

real scalar eresp(real scalar B, real scalar tau, real matrix cf, real scalar bw, real scalar xscale)
{
    real matrix integral, roots, realroots, out

    if (cols(cf) == 1) {
        if (cf[1] <= 0) return(.)
        return(((B * bw) / cf[1]) * xscale)
    }

    integral = polyinteg(cf, 1)
    integral[1] = -polyeval(integral, tau) - B*bw

    roots = polyroots(integral)
    realroots = Re(select(roots, Im(roots) :== 0))
    out = sort(select(realroots, realroots :> tau)', 1)'

    if (cols(out) == 0) return(.)
    else return((out[1] - tau) * xscale)
}

end

