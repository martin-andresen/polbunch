capture program drop polbunchsim
program polbunchsim, rclass
    syntax [, zmin(string) zmax(string) log reps(integer 1) ///
        obs(integer 5000) cutoff(real 1) el(real 0.4) ///
        t0(real 0.2) t1(real 0.6) bw(real 0.01) ///
        bootreps(integer 500) POLynomial(integer 1) ///
        notransform distribution(string) positive ///
        estimator(numlist integer) btype(numlist integer) ///
        clist(string) noisily sample(string)]

    quietly {
        if "`zmin'" == "" local zmin "-."
        if "`zmax'" == "" local zmax "."
        if "`btype'" == "" local btype = 1
        if "`estimator'" == "" local estimator = 3

        if "`clist'" == "" local clist `"noconstant"'
        local numc : word count `clist'
        if `numc' > 2 {
            noi di in red "clist() can contain at most two strings"
            exit 301
        }

        tokenize `clist'
        forvalues i = 1/`numc' {
            if !inlist("``i''", "constant", "noconstant") {
                noi di in red "clist() can contain only constant or noconstant"
                exit 301
            }
        }

        local numb : word count `btype'
        local nume : word count `estimator'
        local numc : word count `clist'
        local numest = `numb' * `nume' * `numc'

        tempname b

        preserve

        capture noisily polbunchgendata z, obs(`obs') cutoff(`cutoff') ///
            el(`el') t0(`t0') t1(`t1') `log' distribution(`distribution')

        local genrc = _rc

        if `genrc' == 0 {
            if "`sample'" != "" drop if !inrange(z, `sample')

            foreach bt of numlist `btype' {
                foreach e of numlist `estimator' {
                    foreach c in `clist' {

                        if "`c'" == "constant" local cval = 1
                        else local cval = 0

                        timer clear
                        timer on 1

                        if `e' == 4 local iff `"if inrange(z, `zmin', `zmax')"'
                        else local iff

                        local rc = 0

                        if `bt' == 0 {
                            capture noisily polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(0) `c' ///
                                `noisily' `notransform' `positive'
                        }
                        else if `bt' == 1 {
                            capture noisily polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(1) `c' ///
                                `noisily' `notransform' `positive'
                        }
                        else if `bt' == 2 {
                            capture noisily bootstrap, reps(`bootreps'): ///
                                polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(0) `c' ///
                                `noisily' `notransform' `positive'
                        }
                        else if `bt' == 3 {
                            capture noisily polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') ///
                                nobayes `c' `noisily' `notransform' `positive'
                        }
                        else if `bt' == 4 {
                            capture noisily polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') ///
                                nozero nobayes `c' `noisily' `notransform' `positive'
                        }
                        else if `bt' == 5 {
                            capture noisily polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') ///
                                `c' `noisily' `notransform' `positive'
                        }
                        else if `bt' == 6 {
                            capture noisily polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') ///
                                nozero `c' `noisily' `notransform' `positive'
                        }
                        else {
                            local rc = 198
                        }

                        if `bt' <= 6 local rc = _rc

                        timer off 1
                        timer list
                        local time = r(t1)
                        timer clear

                        local elast = .
                        local se    = .
                        local p     = .
                        local pmod  = .
                        local novar = 1

                        if `rc' == 0 {
                            capture local elast = _b[bunching:elasticity]
                            if _rc local elast = .

                            capture confirm matrix e(V)
                            if _rc == 0 {
                                capture local se = _se[bunching:elasticity]
                                if _rc == 0 & !missing(`se') & `se' > 0 {
                                    local novar = 0
                                    capture test _b[bunching:elasticity] = `el'
                                    if _rc == 0 local p = r(p)
                                    else local p = .
                                }
                            }

                            capture local pmod = e(p_mod)
                            if _rc local pmod = .
                        }

                        if `numest' > 1 {
                            matrix `b' = nullmat(`b'), `elast', `time'
                            local names `names' e_`bt'_`e'_`cval' t_`bt'_`e'_`cval'

                            if `novar' == 0 {
                                matrix `b' = nullmat(`b'), `se', `p'
                                local names `names' se_`bt'_`e'_`cval' p_`bt'_`e'_`cval'
                            }

                            if !missing(`pmod') {
                                matrix `b' = nullmat(`b'), `pmod'
                                local names `names' p_mod_`bt'_`e'_`cval'
                            }
                        }
                    }
                }
            }
        }

		ereturn clear
       restore

		if `genrc' {
			return scalar failed = 1
			exit
		}

		if `numest' > 1 {
			forvalues j = 1/`=colsof(`b')' {
				local nm : word `j' of `names'
				return scalar `nm' = `b'[1,`j']
			}
		}
		else {
			return scalar e_`btype'_`estimator'_0 = `elast'
			return scalar t_`btype'_`estimator'_0 = `time'
			if `novar' == 0 {
				return scalar se_`btype'_`estimator'_0 = `se'
				return scalar p_`btype'_`estimator'_0  = `p'
			}
			if !missing(`pmod') {
				return scalar p_mod_`btype'_`estimator'_0 = `pmod'
			}
		}
    }
end