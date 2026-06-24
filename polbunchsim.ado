capture program drop polbunchsim
program polbunchsim, eclass
    syntax [, zmin(string) zmax(string) log reps(integer 1) ///
        obs(integer 5000) cutoff(real 1) el(real 0.4) ///
        t0(real 0.2) t1(real 0.6) bw(real 0.01) ///
        bootreps(integer 500) POLynomial(integer 1) ///
        distribution(string) opts(string) ///
        estimator(numlist integer) btype(numlist integer) ///
        clist(string) sample(string)  ///
        est4limits(numlist) limits(numlist)]

    quietly {
        if "`zmin'" == "" local zmin "-."
        if "`zmax'" == "" local zmax "."
        if "`btype'" == "" local btype 1
        if "`estimator'" == "" local estimator 3

        if "`clist'" == "" local clist `"noconstant"'

        local misscode = -1e300

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
        local anyfail = 0

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

                        local modelname m_`bt'_`e'_`cval'_b

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
                                `log' estimator(`e') vce(none) `c' ///
                                 `uselimits' `opts'
                        }
                        else if `bt' == 1 {
                            capture  polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') vce(analytic) `c' ///
                                `opts' `uselimits'  
                        }
                        else if `bt' == 2 {
                            capture bootstrap, reps(`bootreps'): ///
                                polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') vce(none) `c' ///
                                `opts' `uselimits' 
                        }
                        else if `bt' == 3 {
                            capture polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') vce(bootstrap) ///
                                 `c' `opts' `uselimits' 
                        }
                        else if `bt' == 4 {
                            capture polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') vce(bayes) ///
                                 `c' `opts' `uselimits'
                        }
                        else if `bt' == 5 {
                            capture  polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') vce(bootstrap) ///
                                `c' `opts' `uselimits' nozero 
                        }
                        else if `bt' == 6 {
                            capture polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') vce(bayes) ///
                                nozero `c' `opts' `uselimits'
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
                        local novar = 1
                        foreach _tt in wald minimumdistance hausman {
                            local chi2_`_tt' = .
                            local p_`_tt'    = .
                        }

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

                            foreach _tt in wald minimumdistance hausman {
                                capture local chi2_`_tt' = e(chi2_`_tt')
                                if _rc local chi2_`_tt' = .
                                capture local p_`_tt' = e(p_`_tt')
                                if _rc local p_`_tt' = .
                            }

                            if `numest' == 1 {
                                estimates store `esthold'
                            }
                            else {
                                /*
                                    Retain only the coefficients originally posted
                                    by polbunch in the combined e(b). Prefix the
                                    original equation names by the model identifier
                                    so that all coefficient names remain unique.
                                */
                                tempname this_b
                                matrix `this_b' = e(b)

                                local oldnames : colfullnames `this_b'
                                local newnames
                                foreach nm of local oldnames {
                                    gettoken oldeq oldcoef : nm, parse(":")
                                    if "`oldcoef'" == "" {
                                        local oldcoef "`oldeq'"
                                        local oldeq "b"
                                    }
                                    else {
                                        local oldcoef = substr("`oldcoef'", 2, .)
                                    }
                                    local neweq = substr("`modelname'_`oldeq'", 1, 32)
                                    local newnames `newnames' `neweq':`oldcoef'
                                }
                                matrix colnames `this_b' = `newnames'
                                matrix `sim_b' = nullmat(`sim_b'), `this_b'
                            }
                        }

                        /*
                            Quantities computed by polbunchsim, rather than
                            coefficients estimated by polbunch, are returned as
                            model-specific e() scalars. simulate can collect them
                            explicitly as e(<model>_<name>).
                        */
                        local time_post = `time'
                        local p_post    = `p'
                        local se_post   = `se'
                        local rc_post   = `rc'

                        if missing(`time_post') local time_post = `misscode'
                        if missing(`p_post')    local p_post    = `misscode'
                        if missing(`se_post')   local se_post   = `misscode'

                        foreach _tt in wald minimumdistance hausman {
                            local chi2_`_tt'_post = `chi2_`_tt''
                            local p_`_tt'_post    = `p_`_tt''
                            if missing(`chi2_`_tt'_post') local chi2_`_tt'_post = `misscode'
                            if missing(`p_`_tt'_post')    local p_`_tt'_post    = `misscode'
                        }

                        local scalar_names `scalar_names' ///
                            `modelname'_time `modelname'_se `modelname'_p ///
                            `modelname'_chi2_wald `modelname'_p_wald ///
                            `modelname'_chi2_minimumdistance `modelname'_p_minimumdistance ///
                            `modelname'_chi2_hausman `modelname'_p_hausman ///
                            `modelname'_rc

                        local scalar_`modelname'_time             = `time_post'
                        local scalar_`modelname'_se               = `se_post'
                        local scalar_`modelname'_p                = `p_post'
                        local scalar_`modelname'_chi2_wald        = `chi2_wald_post'
                        local scalar_`modelname'_p_wald           = `p_wald_post'
                        local scalar_`modelname'_chi2_minimumdistance = `chi2_minimumdistance_post'
                        local scalar_`modelname'_p_minimumdistance   = `p_minimumdistance_post'
                        local scalar_`modelname'_chi2_hausman     = `chi2_hausman_post'
                        local scalar_`modelname'_p_hausman        = `p_hausman_post'
                        local scalar_`modelname'_rc               = `rc_post'

                        if `rc' local anyfail = 1

                        local final_bt   `bt'
                        local final_e    `e'
                        local final_cval `cval'
                        local final_rc   `rc'
                        local final_time `time'
                    }
                }
            }
        }

        restore

        ereturn clear

        if `genrc' {
            ereturn scalar failed = 1
            ereturn scalar genrc  = `genrc'
            ereturn scalar misscode = `misscode'
            ereturn local cmd "polbunchsim"
            exit
        }

        if `numest' > 1 {
            capture confirm matrix `sim_b'
            if _rc {
                ereturn scalar failed  = 1
                ereturn scalar misscode = `misscode'
                foreach s of local scalar_names {
                    ereturn scalar `s' = `scalar_`s''
                }
                ereturn local cmd "polbunchsim"
                exit
            }

            /*
                e(b) now contains only coefficients that came from the
                constituent polbunch calls. No timing, p-values, or other
                simulation diagnostics are appended to it.
            */
            ereturn post `sim_b'

            foreach s of local scalar_names {
                ereturn scalar `s' = `scalar_`s''
            }

            ereturn scalar failed = `anyfail'
            ereturn scalar misscode = `misscode'
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
                ereturn clear

                foreach s of local scalar_names {
                    ereturn scalar `s' = `scalar_`s''
                }

                ereturn scalar failed  = 1
                ereturn scalar rc      = `final_rc'
                ereturn scalar misscode = `misscode'
                ereturn scalar obs     = `obs'
                ereturn scalar cutoff  = `cutoff'
                ereturn scalar el      = `el'
                ereturn scalar t0      = `t0'
                ereturn scalar t1      = `t1'
                ereturn scalar bw      = `bw'
                ereturn scalar polynomial = `polynomial'
                ereturn scalar bootreps   = `bootreps'

                ereturn local btype "`btype'"
                ereturn local estimator "`estimator'"
                ereturn local clist "`clist'"
                ereturn local cmd "polbunchsim"
                exit
            }

            /*
                Restoring the stored result restores the original polbunch
                e(b), e(V), matrices, macros, and scalars. We then add only
                simulation diagnostics as e() scalars; e(b) is untouched.
            */
            estimates restore `esthold'

            foreach s of local scalar_names {
                ereturn scalar `s' = `scalar_`s''
            }

            ereturn scalar failed = 0
            ereturn scalar misscode = `misscode'

            /*
                Convenient generic aliases for the one-model case.
            */
            ereturn scalar sim_time            = `time_post'
            ereturn scalar sim_se              = `se_post'
            ereturn scalar sim_p               = `p_post'
            ereturn scalar sim_chi2_wald       = `chi2_wald_post'
            ereturn scalar sim_p_wald          = `p_wald_post'
            ereturn scalar sim_chi2_minimumdistance = `chi2_minimumdistance_post'
            ereturn scalar sim_p_minimumdistance    = `p_minimumdistance_post'
            ereturn scalar sim_chi2_hausman    = `chi2_hausman_post'
            ereturn scalar sim_p_hausman       = `p_hausman_post'
            ereturn scalar sim_rc              = `rc_post'

            ereturn scalar obs    = `obs'
            ereturn scalar cutoff = `cutoff'
            ereturn scalar el     = `el'
            ereturn scalar t0     = `t0'
            ereturn scalar t1     = `t1'
            ereturn scalar bw     = `bw'
            ereturn scalar polynomial = `polynomial'
            ereturn scalar bootreps   = `bootreps'

            ereturn local btype "`btype'"
            ereturn scalar estimator=`estimator'
            ereturn local clist "`clist'"
            ereturn local cmd "polbunchsim"
        }
    }
end