	// Helper program for the binned bootstrap
program bssim, eclass 
	syntax, modstr(string) y(varname) estimator(integer) cutoff(real) polynomial(integer) bw(real) h(integer) l(integer)  p(varname) obs(real) [b0(string) estopts(string) t0(numlist min=1 max=1  >=0 <=1) t1(numlist min=1 max=1  >=0 <=1) nl nlcom(string) log constant nopositiveshift boot notransform]
	
	//preserve
	
	replace `y'=.
	loc i=0
	loc factor=0
	while `obs'>0&`i'<`=_N-1' {
		loc ++i
		replace `y'=rbinomial(`obs',`p'/(1-`factor')) in `i'
		loc factor=`factor'+`p'[`i']
		loc obs=`obs'-`y'[`i']
	}
	if `obs'>0 replace `y'=`obs' if `y'==.
	else replace `y'=0 if `y'==.
	
	if "`nl'"=="" reg `y' `modstr', `estopts'
	else nl (`y'=`modstr'), `estopts'
	
	tempname b
	if "`transform'"!="notransform" {
		bunchcalc, estimator(`estimator') polynomial(`polynomial') cutoff(`cutoff') bw(`bw') h(`h') l(`l') b0(`b0') t0(`t0') t1(`t1') `constant' `positiveshift' `log' boot nlcom(`nlcom')
		mat `b'=r(b)
	}
	else mat `b'=e(b)
	ereturn post `b'
	//restore
	//drop `padj'
	end


