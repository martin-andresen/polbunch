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
		/*
			Parse triangular(a,b,c) by locating its parentheses and
			splitting on commas, rather than reading fixed character
			positions (substr(...,12,1) etc.) -- the fixed-position
			version silently misparsed anything but single-digit,
			non-negative a/b/c (e.g. "triangular(-1,10,10)" or
			"triangular(0.5,3,1)" would read nonsense substrings).
			Done outside the "cap" below (with its own clear error
			message) rather than inside it, so a bad argument count
			doesn't just fall through to the generic "use a valid
			distribution" message with no indication of what was
			actually wrong.
		*/
		loc _popen  = strpos("`distribution'","(")
		loc _pclose = strpos("`distribution'",")")
		loc _argstr = substr("`distribution'", `_popen'+1, `_pclose'-`_popen'-1)
		loc _argstr = subinstr("`_argstr'", ","," ",.)
		loc _nargs : word count `_argstr'
		if `_nargs' != 3 {
			noi di as error "triangular(a,b,c) requires exactly 3 numeric arguments; got `_nargs' in distribution(`distribution')."
			exit 198
		}
		loc a : word 1 of `_argstr'
		loc b : word 2 of `_argstr'
		loc c : word 3 of `_argstr'
		cap {
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
	