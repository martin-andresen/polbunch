program polbunchsim, eclass
	syntax , obs(real)  zstar(real) bw(real) pol(real) el(real) t0(real) t1(real) [notransform distribution(string) log noisily nopositiveshift estimator(integer 3)  constant bootreps(integer 0)]
	 {
	polbunchgendata , distribution(`distribution') obs(`obs') zstar(`zstar') el(`el') t0(`t0') t1(`t1') `log'	
	polbunch z, cutoff(`zstar') pol(`pol') bw(`bw') t0(`t0') t1(`t1') `log' `transform' `noisily' `positiveshift' estimator(`estimator') bootreps(`bootreps') `constant'
	tempname b
	mat `b'=e(b)
	ereturn post `b'
	}
end