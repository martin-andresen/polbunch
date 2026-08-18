capture program drop polbunchsim
program polbunchsim, eclass
    syntax [, zmin(string) zmax(string) log reps(integer 1) ///
        obs(integer 5000) cutoff(real 1) el(real 0.4) ///
        t0(real 0.2) t1(real 0.6) bw(real 0.01) ///
        bootreps(integer 500) POLynomial(integer 1) ///
        distribution(string) opts(string) ///
        estimator(numlist integer) btype(numlist integer) ///
        clist(string) sample(string)  ///
        est4limits(numlist) limits(numlist) report(string)]

    quietly {
        if "`zmin'" == "" local zmin "-."
        if "`zmax'" == "" local zmax "."
        if "`btype'" == "" local btype 1
        if "`estimator'" == "" local estimator 3

        if "`clist'" == "" local clist `"noconstant"'

        local misscode = -1e300

        // Split report() into coefficient names (from e(b)) and scalar names
        local _report_coefs ""
        local _report_scals ""
        if "`report'" != "" {
            foreach _item of local report {
                if inlist("`_item'", "hausman", "wald", "minimumdistance") {
                    local _report_scals `_report_scals' chi2_`_item' p_`_item'
                }
                else if inlist("`_item'", ///
                    "chi2_wald","p_wald", ///
                    "chi2_hausman","p_hausman", ///
                    "chi2_minimumdistance","p_minimumdistance", ///
                    "delta_md") {
                    local _report_scals `_report_scals' `_item'
                }
                else {
                    local _report_coefs `_report_coefs' `_item'
                }
            }
        }

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

                        local modelname b`bt'e`e'c`cval'

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

                        local elast    = .
                        local se       = .
                        local p        = .
                        local delta_md = .
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
                            capture local delta_md = e(delta_md)
                            if _rc local delta_md = .

                            if `numest' == 1 {
                                estimates store `esthold'
                            }
                            else {
                                if "`report'" == "" {
                                    /*
                                        No filter: include all e(b) columns with
                                        model-prefixed equation names for uniqueness.
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
                                        /*
                                            simulate builds each result variable's
                                            name from eq+coef together, capped at
                                            32 chars; if we only cap the eq part,
                                            the combined name can still overflow
                                            and simulate silently falls back to
                                            generic names like _sim_43. Reserve
                                            room for the coefficient part first.
                                        */
                                        local _avail = 32 - strlen("`oldcoef'") - 1
                                        if `_avail' < 1 local _avail = 1
                                        local neweq = substr("`modelname'_`oldeq'", 1, `_avail')
                                        local newnames `newnames' `neweq':`oldcoef'
                                    }
                                    matrix colnames `this_b' = `newnames'
                                    matrix `sim_b' = nullmat(`sim_b'), `this_b'
                                }
                                else {
                                    /*
                                        report() specified: flat naming {item}_{bt}_{e}_{cval}
                                        with empty equations, combining coef and scalar columns.
                                    */
                                    tempname this_b

                                    if "`_report_coefs'" != "" {
                                        tempname _cof_b
                                        mata: _pbsim_filter_b("e(b)", "`_cof_b'", "`_report_coefs'")
                                        capture confirm matrix `_cof_b'
                                        if !_rc {
                                            local _nc = colsof(`_cof_b')
                                            local _news ""
                                            forvalues _j = 1/`_nc' {
                                                local _allnm : colfullnames `_cof_b'
                                                local _nm : word `_j' of `_allnm'
                                                gettoken _eq _cn : _nm, parse(":")
                                                if "`_cn'" == "" local _cn "`_eq'"
                                                else local _cn = substr("`_cn'", 2, .)
                                                local _suffix _`bt'_`e'_`cval'
                                                local _avail = 32 - strlen("`_suffix'")
                                                if `_avail' < 1 local _avail = 1
                                                local _news `_news' `=substr("`_cn'", 1, `_avail')'`_suffix'
                                            }
                                            mata: _pbsim_set_flat_stripe("`_cof_b'", "`_news'")
                                            matrix `this_b' = `_cof_b'
                                        }
                                    }

                                    foreach _sc of local _report_scals {
                                        local _scval = ``_sc''
                                        if missing(`_scval') continue
                                        tempname _scol
                                        matrix `_scol' = (`_scval')
                                        mata: _pbsim_set_flat_stripe("`_scol'", "`_sc'_`bt'_`e'_`cval'")
                                        capture confirm matrix `this_b'
                                        if _rc matrix `this_b' = `_scol'
                                        else  matrix `this_b' = `this_b', `_scol'
                                    }

                                    capture confirm matrix `this_b'
                                    if !_rc matrix `sim_b' = nullmat(`sim_b'), `this_b'
                                }
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

                        local `modelname'_time                 = `time_post'
                        local `modelname'_se                   = `se_post'
                        local `modelname'_p                    = `p_post'
                        local `modelname'_chi2_wald            = `chi2_wald_post'
                        local `modelname'_p_wald               = `p_wald_post'
                        local `modelname'_chi2_minimumdistance = `chi2_minimumdistance_post'
                        local `modelname'_p_minimumdistance    = `p_minimumdistance_post'
                        local `modelname'_chi2_hausman         = `chi2_hausman_post'
                        local `modelname'_p_hausman            = `p_hausman_post'
                        local `modelname'_rc                   = `rc_post'

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
                    ereturn scalar `s' = ``s''
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
                ereturn scalar `s' = ``s''
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
                    ereturn scalar `s' = ``s''
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
                simulation diagnostics as e() scalars.
            */
            estimates restore `esthold'

            if "`report'" != "" {
                tempname combined_b

                // Coefficient columns: preserve polbunch equation names
                if "`_report_coefs'" != "" {
                    tempname _filt_b
                    mata: _pbsim_filter_b("e(b)", "`_filt_b'", "`_report_coefs'")
                    capture confirm matrix `_filt_b'
                    if !_rc {
                        capture confirm matrix `combined_b'
                        if _rc matrix `combined_b' = `_filt_b'
                        else  matrix `combined_b' = `combined_b', `_filt_b'
                    }
                }

                // Scalar columns: eq = test name, coef = stat type
                foreach _sc of local _report_scals {
                    local _scval = ``_sc''
                    if missing(`_scval') continue
                    local _sc_eq "test"
                    local _sc_cn "`_sc'"
                    foreach _test in wald hausman minimumdistance {
                        if "`_sc'" == "chi2_`_test'" {
                            local _sc_eq `_test'
                            local _sc_cn chi2
                        }
                        if "`_sc'" == "p_`_test'" {
                            local _sc_eq `_test'
                            local _sc_cn p
                        }
                    }
                    if "`_sc'" == "delta_md" {
                        local _sc_eq minimumdistance
                        local _sc_cn delta
                    }
                    tempname _scol
                    matrix `_scol' = (`_scval')
                    matrix colnames `_scol' = `_sc_cn'
                    matrix coleq   `_scol' = `_sc_eq'
                    capture confirm matrix `combined_b'
                    if _rc matrix `combined_b' = `_scol'
                    else  matrix `combined_b' = `combined_b', `_scol'
                }

                // Post combined_b; include V only when there are no scalar columns
                capture confirm matrix `combined_b'
                if !_rc {
                    if "`_report_scals'" == "" {
                        capture confirm matrix e(V)
                        if !_rc {
                            tempname _filt_V
                            mata: _pbsim_filter_V("e(V)", "`_filt_V'", "`_report_coefs'")
                            capture confirm matrix `_filt_V'
                            if !_rc ereturn post `combined_b' `_filt_V'
                            else    ereturn post `combined_b'
                        }
                        else ereturn post `combined_b'
                    }
                    else ereturn post `combined_b'
                }
            }

            foreach s of local scalar_names {
                ereturn scalar `s' = ``s''
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

mata:
void _pbsim_filter_b(string scalar src, string scalar dst, string scalar keep_str) {
    real matrix B
    string matrix stripe
    real rowvector idx
    real scalar j
    string colvector keep

    B      = st_matrix(src)
    stripe = st_matrixcolstripe(src)
    keep   = tokens(keep_str)'

    idx = J(1, 0, .)
    for (j = 1; j <= rows(stripe); j++) {
        if (anyof(keep, stripe[j, 2])) idx = (idx, j)
    }
    if (length(idx) == 0) return

    st_matrix(dst, B[., idx])
    st_matrixcolstripe(dst, stripe[idx',])
    st_matrixrowstripe(dst, st_matrixrowstripe(src))
}

void _pbsim_filter_V(string scalar src, string scalar dst, string scalar keep_str) {
    real matrix V
    string matrix stripe
    real rowvector idx
    real scalar j
    string colvector keep

    V      = st_matrix(src)
    stripe = st_matrixcolstripe(src)
    keep   = tokens(keep_str)'

    idx = J(1, 0, .)
    for (j = 1; j <= rows(stripe); j++) {
        if (anyof(keep, stripe[j, 2])) idx = (idx, j)
    }
    if (length(idx) == 0) return

    st_matrix(dst, V[idx, idx])
    st_matrixcolstripe(dst, stripe[idx',])
    st_matrixrowstripe(dst, stripe[idx',])
}

void _pbsim_set_flat_stripe(string scalar mat, string scalar names_str) {
    string colvector names
    names = tokens(names_str)'
    st_matrixcolstripe(mat, (J(rows(names), 1, ""), names))
}
end