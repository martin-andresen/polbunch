capture program drop polbunchsim
program polbunchsim, eclass
    syntax [, zmin(string) zmax(string) log reps(integer 1) ///
        obs(integer 5000) cutoff(real 1) el(real 0.4) ///
        t0(real 0.2) t1(real 0.6) bw(real 0.01) ///
        bootreps(integer 500) POLynomial(integer 1) ///
        notransform distribution(string) positive ///
        estimator(numlist integer) btype(numlist integer) ///
        clist(string) noisily sample(string) limits(numlist) ///
        est4limits(numlist) minimumdistance norankred]

    quietly {
        if "`zmin'" == "" local zmin "-."
        if "`zmax'" == "" local zmax "."
        if "`btype'" == "" local btype 1
        if "`estimator'" == "" local estimator 3

        if "`clist'" == "" local clist `"noconstant"'

        local numc : word count `clist'
        if `numc' > 2 {
            noi di as error "clist() can contain at most two strings"
            exit 301
        }

        tokenize `clist'
        forvalues i = 1/`numc' {
            if !inlist("``i''", "constant", "noconstant") {
                noi di as error "clist() can contain only constant or noconstant"
                exit 301
            }
        }

        local numb : word count `btype'
        local nume : word count `estimator'
        local numc : word count `clist'
        local numest = `numb' * `nume' * `numc'

        tempname sim_b sim_extra sim_oldb sim_oldV sim_newb sim_newV
        tempname esthold

        preserve

        capture noisily polbunchgendata z, obs(`obs') cutoff(`cutoff') ///
            el(`el') t0(`t0') t1(`t1') `log' distribution(`distribution')

        local genrc = _rc

        if `genrc' == 0 {
            if "`sample'" != "" {
                drop if !inrange(z, `sample')
            }

            foreach bt of numlist `btype' {
                foreach e of numlist `estimator' {
                    foreach c in `clist' {

                        if "`c'" == "constant" local cval = 1
                        else local cval = 0

                        local modelname m_`bt'_`e'_`cval'

                        timer clear
                        timer on 1

                        if `e' == 4 local iff `"if inrange(z, `zmin', `zmax')"'
                        else local iff

                        if `e' == 4 & "`est4limits'" != "" {
                            local uselimits limits(`est4limits')
                        }
                        else {
                            local uselimits limits(`limits')
                        }

                        local rc = 0

                        if `bt' == 0 {
                            capture polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(0) `c' ///
                                `noisily' `notransform' `positive' ///
                                `uselimits' `minimumdistance' `rankred'
                        }
                        else if `bt' == 1 {
                            capture polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(1) `c' ///
                                `noisily' `notransform' `positive' ///
                                `uselimits' `minimumdistance' `rankred'
                        }
                        else if `bt' == 2 {
                            capture bootstrap, reps(`bootreps'): ///
                                polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(0) `c' ///
                                `noisily' `notransform' `positive' ///
                                `uselimits' `minimumdistance' `rankred'
                        }
                        else if `bt' == 3 {
                            capture polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') ///
                                nobayes `c' `noisily' `notransform' `positive' ///
                                `uselimits' `minimumdistance' `rankred'
                        }
                        else if `bt' == 4 {
                            capture polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') ///
                                nozero nobayes `c' `noisily' `notransform' ///
                                `positive' `uselimits' `minimumdistance' `rankred'
                        }
                        else if `bt' == 5 {
                            capture  polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') ///
                                `c' `noisily' `notransform' `positive' ///
                                `uselimits' `minimumdistance' `rankred'
                        }
                        else if `bt' == 6 {
                            capture polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') ///
                                nozero `c' `noisily' `notransform' `positive' ///
                                `uselimits' `minimumdistance' `rankred'
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

                            if `numest' == 1 {
                                estimates store `esthold'
                            }
                        }

                        if `numest' > 1 {
                            matrix `sim_extra' = (`elast', `time')
                            matrix colnames `sim_extra' = elasticity time
                            matrix coleq `sim_extra' = `modelname' `modelname'

                            matrix `sim_b' = nullmat(`sim_b'), `sim_extra'

                            if `novar' == 0 {
                                matrix `sim_extra' = (`se', `p')
                                matrix colnames `sim_extra' = se p
                                matrix coleq `sim_extra' = `modelname' `modelname'

                                matrix `sim_b' = nullmat(`sim_b'), `sim_extra'
                            }

                            if !missing(`pmod') {
                                matrix `sim_extra' = (`pmod')
                                matrix colnames `sim_extra' = p_mod
                                matrix coleq `sim_extra' = `modelname'

                                matrix `sim_b' = nullmat(`sim_b'), `sim_extra'
                            }
                        }

                        local final_bt   `bt'
                        local final_e    `e'
                        local final_cval `cval'
                        local final_rc   `rc'
                    }
                }
            }
        }

        restore

        ereturn clear

        if `genrc' {
            ereturn scalar failed = 1
            ereturn scalar genrc  = `genrc'
            ereturn local cmd "polbunchsim"
            exit
        }

        if `numest' > 1 {
            capture confirm matrix `sim_b'
            if _rc {
                ereturn scalar failed = 1
                ereturn local cmd "polbunchsim"
                exit
            }

            ereturn post `sim_b'

            ereturn scalar failed = 0
            ereturn scalar numest = `numest'
            ereturn scalar obs    = `obs'
            ereturn scalar cutoff = `cutoff'
            ereturn scalar el     = `el'
            ereturn scalar t0     = `t0'
            ereturn scalar t1     = `t1'
            ereturn scalar bw     = `bw'
            ereturn scalar polynomial = `polynomial'
            ereturn scalar bootreps   = `bootreps'

            ereturn local btype "`btype'"
            ereturn local estimator "`estimator'"
            ereturn local clist "`clist'"
            ereturn local cmd "polbunchsim"
        }
        else {
            if `final_rc' {
                matrix `sim_b' = (., `time')
                matrix colnames `sim_b' = m_`final_bt'_`final_e'_`final_cval':elasticity ///
                                          m_`final_bt'_`final_e'_`final_cval':time

                ereturn post `sim_b'

                ereturn scalar failed = 1
                ereturn scalar rc     = `final_rc'
                ereturn scalar time   = `time'
                ereturn local cmd "polbunchsim"
                exit
            }

            estimates restore `esthold'

            local escalars  : e(scalars)
            local emacros   : e(macros)
            local ematrices : e(matrices)

            foreach s of local escalars {
                local hold_scalar_`s' = e(`s')
            }

            foreach m of local emacros {
                local hold_macro_`m' `"`e(`m')'"'
            }

            local mh 0
            foreach M of local ematrices {
                if !inlist("`M'", "b", "V") {
                    local ++mh
                    tempname hmat`mh'
                    matrix `hmat`mh'' = e(`M')
                    local hmat_name_`mh' "`M'"
                }
            }
            local hmat_count `mh'

            matrix `sim_oldb' = e(b)

            capture confirm matrix e(V)
            local hasV = (_rc == 0)

            if `hasV' {
                matrix `sim_oldV' = e(V)
            }

            matrix `sim_extra' = (`time')
            matrix colnames `sim_extra' = polbunchsim:time

            if `novar' == 0 {
                matrix `sim_extra' = `sim_extra', (`p')
                matrix colnames `sim_extra' = polbunchsim:time polbunchsim:p
            }

            if !missing(`pmod') {
                matrix `sim_extra' = `sim_extra', (`pmod')

                if `novar' == 0 {
                    matrix colnames `sim_extra' = polbunchsim:time polbunchsim:p polbunchsim:p_mod
                }
                else {
                    matrix colnames `sim_extra' = polbunchsim:time polbunchsim:p_mod
                }
            }

            matrix `sim_newb' = `sim_oldb', `sim_extra'

            if `hasV' {
                local k_old = colsof(`sim_oldb')
                local k_new = colsof(`sim_newb')

                matrix `sim_newV' = J(`k_new', `k_new', 0)

                forvalues ii = 1/`k_old' {
                    forvalues jj = 1/`k_old' {
                        matrix `sim_newV'[`ii', `jj'] = `sim_oldV'[`ii', `jj']
                    }
                }
				
				loc names: colfullnames `sim_newb'

                matrix colnames `sim_newV' = `names'
                matrix rownames `sim_newV' = `names'
				
                ereturn post `sim_newb' `sim_newV'
            }
            else {
                ereturn post `sim_newb'
            }

            foreach s of local escalars {
                if !inlist("`s'", "rank") {
                    capture ereturn scalar `s' = `hold_scalar_`s''
                }
            }

            foreach m of local emacros {
                capture ereturn local `m' `"`hold_macro_`m''"'
            }

            forvalues hh = 1/`hmat_count' {
                local M "`hmat_name_`hh''"
                capture ereturn matrix `M' = `hmat`hh''
            }

            ereturn scalar failed     = 0
            ereturn scalar sim_time   = `time'
            ereturn scalar sim_p      = `p'
            ereturn scalar sim_p_mod  = `pmod'

            ereturn scalar obs    = `obs'
            ereturn scalar cutoff = `cutoff'
            ereturn scalar el     = `el'
            ereturn scalar t0     = `t0'
            ereturn scalar t1     = `t1'
            ereturn scalar bw     = `bw'
            ereturn scalar polynomial = `polynomial'
            ereturn scalar bootreps   = `bootreps'

            ereturn local btype "`btype'"
            ereturn local estimator "`estimator'"
            ereturn local clist "`clist'"
            ereturn local cmd "polbunchsim"
			
        }
		
    }
	ereturn display
end