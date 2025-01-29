program polbunchsim, eclass
	syntax newvarname, obs(real) cutoff(real) bw(real) pol(real) el(real) t0(real) t1(real) [notransform distribution(string) log noisily nopositiveshift estimator(integer 3) constant bootreps(integer 0) notest nobayes nozero]
	 {

	polbunchgendata `varlist', distribution(`distribution') obs(`obs') cutoff(`cutoff') el(`el') t0(`t0') t1(`t1') `log'	
	timer on 99
	if `bootreps'>0 polbunch `varlist', cutoff(`cutoff') pol(`pol') bw(`bw') t0(`t0') t1(`t1') `log' `transform' `noisily' `positiveshift' estimator(`estimator') bootreps(`bootreps') `constant' `test' `bayes' `zero'
	else if `bootreps'==0  polbunch `varlist', cutoff(`cutoff') pol(`pol') bw(`bw') t0(`t0') t1(`t1') `log' `transform' `noisily' `positiveshift' estimator(`estimator') bootreps(0) `constant'
	else bootstrap, reps(`=-`bootreps''): polbunch `varlist', cutoff(`cutoff') pol(`pol') bw(`bw') t0(`t0') t1(`t1') `log' `transform' `noisily' `positiveshift' estimator(`estimator') bootreps(0) `constant'
	if `bootreps'>=1 loc p_mod=e(p_mod)
	else loc p_mod=.
	tempname b V
	mat `b'=e(b)
	if `bootreps'==0 ereturn post `b'
	else {
		mat `V'=e(V)
		ereturn post `b' `V'
	}
	timer off 99
	timer list
	timer clear
	ereturn scalar time=`r(t99)'
	if `bootreps'!=0 {
		test elasticity=`el'
		ereturn scalar p=r(p)
		}
	ereturn scalar p_mod=`p_mod'
	ereturn local cmd "polbunch"
	ereturn scalar polynomial=e(polynomial)
	drop `varlist'
	}
	
end