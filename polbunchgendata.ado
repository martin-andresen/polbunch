program polbunchgendata, rclass
	syntax newvarname [,obs(integer 5000) cutoff(real 1) el(real 0.4) t0(real 0.2) t1(real 0.6) distribution(string) log buncherror(string)]
	
	quietly {
	clear
	set obs `obs'
	
	if "`buncherror'" !="" {
		tmp error
		cap gen `error'=`buncherror'
		if _rc!=0 {
			noi di as error "Option buncherror must be an expression that can generate a Stata variable - it should indicate the optimization friction or error for bunchers. Try e.g. +rnormal() for an additive normally distributed error or *exp(rnormal()) for a multiplicative normal error. "
			exit 301
		}
	}
	if !inrange(`t0',0,1)|!inrange(`t1',0,1)|`t1'<`t0'|`el'<0 {
		noi di as error "tax rates t1 and t0 must be between 0 and 1, t1>t0, and the elasticity must be nonnegative. "
	}
	tempname r
	if "`distribution'"=="" loc distribution triangular(0,3,0)
	if strpos("`distribution'","triangular")>0 {
		cap {
		loc a=substr("`distribution'",12,1)
		loc b=substr("`distribution'",14,1)
		loc c=substr("`distribution'",16,1)
		gen double `r'=runiform()
		gen double `varlist'=`a'+sqrt(`r'*(`b'-`a')*(`c'-`a')) if `r'<(`c'-`a')/(`b'-`a')
		replace `varlist'=`b'- sqrt((1-`r')*(`b'-`a')*(`b'-`c')) if `r'>(`c'-`a')/(`b'-`a')
		}
	}
	else cap gen double `varlist'=`distribution'
	if _rc!=0 {
		noi di as error "Use a valid stata random number distribution function (with parameters) in distribution(). Additionally, triangular(a,b,c) is allowed."
	}
	
	su `varlist'
	if !inrange(`cutoff',r(min),r(max)) {
		noi di as error "Cutoff=`cutoff' is not within the support specified by the distribution in distribution()."
		exit 301
	}
	if "`log'"=="log" {
		replace `varlist'=`cutoff'`buncherror' if inrange(`varlist',`cutoff',`cutoff'+`el'*(ln(1-`t0')-ln(1-`t1')))
		replace `varlist'=`varlist'-`el'*(ln(1-`t0')-ln(1-`t1')) if `varlist'>`cutoff'
	}
	else {
		replace `varlist'=`cutoff'`buncherror' if inrange(`varlist',`cutoff',`cutoff'*((1-`t0')/(1-`t1'))^`el')
		replace `varlist'=`varlist'*((1-`t1')/(1-`t0'))^`el' if `varlist'>`cutoff'
	}
	}
	
end
	