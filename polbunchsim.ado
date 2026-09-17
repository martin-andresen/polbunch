capture program drop polbunchsim
program polbunchsim, eclass
    syntax [, zmin(string) zmax(string) log reps(integer 1) ///
        obs(integer 5000) cutoff(real 1) ELasticity(string) ///
        t0(real 0.2) t1(real 0.6) bw(real 0.01) ///
        INCOMEeffect(string) buncherror(string) ///
        bootreps(integer 500) POLynomial(integer 1) ///
        distribution(string) opts(string) ///
        estimator(numlist integer) btype(numlist integer) ///
        clist(string) sample(string)  ///
        est4limits(numlist) limits(numlist) report(string) SCALARsonly DEBUG ///
        PANel(numlist min=3 max=4) B0(real -999) E0(real -999) ///
        HEAP(numlist min=2 max=3) INATTention(real 0) WELfriction(real 0) ///
        PERMStep(real 0) PERMMaxcutoffs(integer 60)]

    /*
        panel(n T rho)  -- pooled-panel DGP, passed straight to
        polbunchgendata (obs is overridden to n*T).  btype 8 (vce(cluster
        pid)) and btype 10 (bootstrap, cluster(pid) prefix) require it.

        b0(#) / e0(#) -- the population excess-mass count and true
        elasticity, for the coverage / size tests.  Compute once from a
        large polbunchgendata draw and pass them in; otherwise b0
        falls back to the realised bunch count of each replication (a
        noisy target -- fine for a smoke test, not for coverage numbers)
        and e0 falls back to the homogeneous elasticity() value.  e0(0) is a
        valid value (size cell); the sentinel for "not given" is < 0.

        permstep(#) / permmaxcutoffs(#) -- passed straight to
        polbunch_permute's step()/maxcutoffs() for btype 12 (see below).
        permstep(0) (the default) resolves to polbunch_permute's own
        default, the fitted model's bw.
    */

    /*
        scalarsonly : single-estimator runs only.  Skip -estimates
        restore- and post NOTHING to e(b); return a fixed set of
        e(sim_*) scalars (genuine missing when a quantity was not
        computed) that is byte-identical in structure across every
        replication -- so -simulate- / -parallel sim- never trip over
        a replication whose underlying polbunch fit failed.
    */

    /*
        elasticity() takes EITHER a nonnegative number or a Stata expression
        (drawn per observation), passed straight to polbunchgendata --
        as do incomeeffect() and buncherror().  When elasticity() is an
        expression there is no single "true" elasticity: the realised
        mean r(el_mean) from polbunchgendata is used as the target of
        the elasticity coverage test and returned in e(eltrue).
    */

    /*
        debug: the internal polbunch calls below are wrapped in a bare
        "capture" (not "capture noisily"), so a failed replication's real
        Stata/Mata error text is normally unrecoverable from outside --
        even calling polbunchsim itself noisily reveals nothing, since
        the inner capture swallows it regardless. With debug specified,
        those internal calls run noisily (still capture'd, so a single
        failed replication still doesn't abort the caller), so the
        actual error surfaces in the log instead of only a bare rc.
    */
    local dbgpfx = cond("`debug'" != "", "noisily", "")

    quietly {
        if "`zmin'" == "" local zmin "-."
        if "`zmax'" == "" local zmax "."
        if "`btype'" == "" local btype 1
        if "`estimator'" == "" local estimator 3
        if "`elasticity'" == "" local elasticity 0.4
        local elnum = real("`elasticity'")

        if "`clist'" == "" local clist `"noconstant"'

        local misscode = -1e300

        // panel: obs is n*T
        local panelopt ""
        if "`panel'" != "" {
            tokenize `panel'
            local obs = `1' * `2'
            local panelopt panel(`panel')
        }

        // Split report() into coefficient names (from e(b)) and scalar names
        local _report_coefs ""
        local _report_scals ""
        /*
            inlist() for strings caps out at 10 comparison values -- this
            list has 11, so membership is checked with strpos() against a
            padded, space-delimited list instead.
        */
        local _scalarnames " chi2_wald p_wald chi2_hausman p_hausman chi2_minimumdistance p_minimumdistance delta_md time se p rc "
        if "`report'" != "" {
            foreach _item of local report {
                if inlist("`_item'", "hausman", "wald", "minimumdistance") {
                    local _report_scals `_report_scals' chi2_`_item' p_`_item'
                }
                else if strpos("`_scalarnames'", " `_item' ") > 0 {
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

        /*
            Only stamp a dimension into the model tag when it actually
            varies in this call, so e.g. a run that only sweeps
            estimator()/clist() gets tags like "e3c0" instead of the
            noisier "b1e3c0" -- and a run that only sweeps btype() gets
            "b1", "b2", etc.
        */
        local vary_bt = (`numb' > 1)
        local vary_e  = (`nume' > 1)
        local vary_c  = (`numc' > 1)

        tempname sim_b sim_V sim_extra sim_oldb sim_oldV sim_newb sim_newV
        tempname esthold

        preserve
        local anyfail = 0

        /*
            polbunchgendata's `syntax newvarname' confirms z does not
            already exist BEFORE any of its body (including its own
            -clear-) runs -- so if the caller's dataset happens to have a
            variable named z (e.g. left behind by a prior polbunchgendata/
            polbunch call, as under -simulate-, which preserves whatever
            was in memory at the moment it was invoked as the per-
            replication baseline), every replication fails with "variable
            z already defined" before generating anything. Safe
            unconditionally: we're inside -preserve- and about to -clear-
            and regenerate everything anyway.
        */
        clear

        local heapopt ""
        if "`heap'" != "" local heapopt heap(`heap')

        capture noisily polbunchgendata z, obs(`obs') cutoff(`cutoff') ///
            elasticity(`elasticity') t0(`t0') t1(`t1') `log' distribution(`distribution') ///
            incomeeffect(`incomeeffect') buncherror(`buncherror') `panelopt' `heapopt' ///
            inattention(`inattention') welfriction(`welfriction')

        local genrc = _rc
        local btrue_realized = .
        if `genrc' == 0 {
            local eltrue = cond(`e0' >= 0, `e0', r(el_mean))
            local ietrue = r(incomeeffect)
            local btrue_realized = r(n_bunchers)
        }
        if "`eltrue'" == "" local eltrue = cond(`e0' >= 0, `e0', `elnum')
        if "`eltrue'" == "" local eltrue = .
        if "`ietrue'" == "" local ietrue = .

        // population excess-mass count for the coverage test
        local btrue_use = cond(`b0' >= 0, `b0', `btrue_realized')

        if `genrc' == 0 {
            if "`sample'" != "" {
                drop if !inrange(z, `sample')
            }

            foreach bt of numlist `btype' {
                foreach e of numlist `estimator' {
                    foreach c in `clist' {

                        if "`c'" == "constant" local cval = 1
                        else local cval = 0

                        local modelname ""
                        if `vary_bt' local modelname "`modelname'b`bt'"
                        if `vary_e'  local modelname "`modelname'e`e'"
                        if `vary_c'  local modelname "`modelname'c`cval'"

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
                        local sim_gof_p        = .
                        local sim_gof_deviance = .
                        local sim_gof_p_below  = .
                        local sim_polyused     = .
                        local permrc           = .
                        local permp            = .
                        local permobserved     = .
                        local permobservedse   = .
                        local permnplacebo     = .
                        local permntried       = .

                        if `bt' == 0 {
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') vce(none) `c' ///
                                 `uselimits' `opts'
                        }
                        else if `bt' == 1 {
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') vce(analytic) `c' ///
                                `opts' `uselimits'
                            /*
                                polbunch silently LOWERS the polynomial order
                                (by 1, possibly repeatedly) when the two
                                one-sided polynomials are rank-deficient at
                                the requested order (its own "Note: Polynomial
                                order lowered..." -- see polbunch.ado). e(),
                                not the requested `polynomial' local, is the
                                only place the ACTUAL fitted degree survives,
                                so it's captured here every rep -- otherwise a
                                sweep over polynomial() (e.g. sim_assumption_
                                tests.do's Panel B) can silently re-fit a
                                lower K than requested with no way to tell.
                            */
                            if _rc == 0 {
                                capture confirm scalar e(polynomial)
                                if !_rc local sim_polyused = e(polynomial)
                            }
                            /*
                                Reference-region GoF/deviance is captured for
                                EVERY estimator, not just 0. For estimator 0
                                (unrestricted: h0, h1 fit separately, no
                                nested restriction to run a Wald/Hausman/MD
                                test against) this is the ONLY assumption
                                test available, and A3 alone at that (see
                                test.tex's "below only" paragraph). For
                                estimators 1/2/3, this is POOLED deviance on
                                the RESTRICTED fit -- the same overall null as
                                the omnibus Wald/MD test (that estimator's
                                cross-kink restriction, given A2/A3), just
                                viewed through the many-bin reference-region
                                fit rather than the few-parameter structural
                                contrast; the two are complementary (weak vs.
                                strong against smooth vs. localized
                                departures respectively -- confirmed by
                                direct simulation, sim_assumption_tests.do's
                                panels B/C). Captured here, immediately after
                                the fit, so no later diagnostic call
                                (Hausman, bootstrap, ...) can overwrite e()
                                first.
                            */
                            if _rc == 0 {
                                capture estat gof
                                if _rc == 0 {
                                    local sim_gof_p        = e(deviance_p)
                                    local sim_gof_deviance = e(deviance)
                                    /*
                                        below-cutoff-only deviance p-value: a
                                        PURE A3 test (h0's own fit against
                                        clean reference data, untouched by any
                                        behavioral response) -- unlike the
                                        pooled deviance_p above, which mixes
                                        in h1's fit (A3 + A1/A2). Exposed
                                        separately since a DGP that only
                                        violates A3 (e.g. sim_assumption_
                                        tests.do's panel B) should show more
                                        power here than in the pooled version.
                                    */
                                    capture confirm scalar e(deviance_below_p)
                                    if !_rc local sim_gof_p_below = e(deviance_below_p)
                                }
                            }
                        }
                        else if `bt' == 2 {
                            /* individual (person-year) nonparametric bootstrap,
                               Saez-2010 style: resample rows, re-bin, re-fit */
                            capture `dbgpfx' bootstrap ///
                                _bs_el = _b[bunching:elasticity] ///
                                _bs_nb = _b[bunching:number_bunchers], ///
                                reps(`bootreps') nodots: ///
                                polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') vce(none) `c' ///
                                `opts' `uselimits'
                        }
                        else if `bt' == 3 {
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') vce(bootstrap) ///
                                 `c' `opts' `uselimits'
                        }
                        else if `bt' == 4 {
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') vce(bayes) ///
                                 `c' `opts' `uselimits'
                        }
                        else if `bt' == 5 {
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') vce(bootstrap) ///
                                `c' `opts' `uselimits' nozero
                        }
                        else if `bt' == 6 {
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') vce(bayes) ///
                                nozero `c' `opts' `uselimits'
                        }
                        else if `bt' == 7 {
                            /* analytic Eicker-White sandwich */
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') vce(robust) `c' ///
                                `opts' `uselimits'
                        }
                        else if `bt' == 8 {
                            /* analytic cluster-robust (needs panel -> pid) */
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') vce(cluster pid) `c' ///
                                `opts' `uselimits'
                        }
                        else if `bt' == 9 {
                            /* Chetty/CFOP residual bootstrap over bins */
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') ///
                                vce(bootstrap, residual) `c' `opts' `uselimits'
                        }
                        else if `bt' == 10 {
                            /* person-clustered nonparametric bootstrap (benchmark
                               for vce(cluster) under pooled-panel dependence) */
                            capture `dbgpfx' bootstrap ///
                                _bs_el = _b[bunching:elasticity] ///
                                _bs_nb = _b[bunching:number_bunchers], ///
                                reps(`bootreps') cluster(pid) nodots: ///
                                polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') vce(none) `c' ///
                                `opts' `uselimits'
                        }
                        else if `bt' == 11 {
                            /* quasi-multinomial: conventional V scaled by the
                               Pearson dispersion phi-hat (Dirichlet-multinomial
                               / Kish design-effect fallback when M is
                               unavailable -- histogram-only) */
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') vce(conventional) scale(x2) `c' ///
                                `opts' `uselimits'
                        }
                        else if `bt' == 12 {
                            /* polbunch_permute: placebo-cutoff permutation
                               inference (r(p) against H0: no bunching
                               anywhere), layered on top of an ordinary
                               vce(analytic) fit so sim_b/sim_se_b/sim_el
                               stay comparable to btype 1 in the same cell.
                               polbunch_permute performs its own separate,
                               always-allownegative internal refit at the
                               true cutoff -- intentionally decoupled from
                               this headline point estimate -- so its
                               output is captured into sim_perm_* below,
                               not folded into sim_b/sim_se_b. */
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`polynomial') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') vce(analytic) `c' ///
                                `opts' `uselimits'
                            /* the base fit's own rc drives elast/b_bunch
                               extraction below (via `rc', set from this
                               local right after the branch) -- kept
                               separate from polbunch_permute's own rc
                               (permrc) so a permute failure never discards
                               an otherwise-successful headline point
                               estimate. */
                            local basefit_rc = _rc
                            if `basefit_rc' == 0 {
                                /* polbunch_permute refits internally (the
                                   verification/observed refit, then every
                                   placebo); protect this headline fit's
                                   e() from whatever it leaves active if it
                                   exits mid-way, so elast/b_bunch below
                                   always come from THIS fit regardless of
                                   how the permute call itself fares. */
                                tempname _basehold
                                quietly estimates store `_basehold'
                                capture `dbgpfx' polbunch_permute, ///
                                    target(bunching:excess_mass) ///
                                    step(`permstep') maxcutoffs(`permmaxcutoffs') nodots
                                if _rc == 0 {
                                    local permrc          = 0
                                    local permp            = r(p)
                                    local permobserved     = r(observed)
                                    local permobservedse   = r(observed_se)
                                    local permnplacebo     = r(n_placebo)
                                    local permntried       = r(n_placebo_tried)
                                }
                                else local permrc = _rc
                                quietly estimates restore `_basehold'
                                estimates drop `_basehold'
                            }
                        }
                        else {
                            local rc = 198
                        }

                        if `bt' == 12 local rc = `basefit_rc'
                        else if `bt' <= 11 local rc = _rc

                        timer off 1
                        timer list
                        local time = r(t1)
                        timer clear

                        local elast    = .
                        local se       = .
                        local p        = .
                        local b_bunch  = .
                        local se_bunch = .
                        local p_bunch  = .
                        local delta_md = .
                        local novar = 1
                        foreach _tt in wald minimumdistance hausman {
                            local chi2_`_tt' = .
                            local p_`_tt'    = .
                        }

                        /* coefficient names: the bootstrap-prefix btypes
                           (2, 10) collect the two targets as _bs_el / _bs_nb;
                           everything else keeps polbunch's own equation names */
                        if inlist(`bt', 2, 10) {
                            local nm_el "_bs_el"
                            local nm_nb "_bs_nb"
                        }
                        else {
                            local nm_el "bunching:elasticity"
                            local nm_nb "bunching:number_bunchers"
                        }

                        if `rc' == 0 {
                            capture local elast = _b[`nm_el']
                            if _rc local elast = .
                            capture local b_bunch = _b[`nm_nb']
                            if _rc local b_bunch = .

                            capture confirm matrix e(V)
                            if _rc == 0 {
                                capture local se = _se[`nm_el']
                                if _rc == 0 & !missing(`se') & `se' > 0 {
                                    local novar = 0
                                    capture test _b[`nm_el'] = `eltrue'
                                    if _rc == 0 local p = r(p)
                                    else local p = .
                                }
                                capture local se_bunch = _se[`nm_nb']
                                if _rc == 0 & !missing(`se_bunch') & `se_bunch' > 0 & `btrue_use' < . {
                                    capture test _b[`nm_nb'] = `btrue_use'
                                    if _rc == 0 local p_bunch = r(p)
                                    else local p_bunch = .
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
                        }

                        /*
                            Quantities computed by polbunchsim, rather than
                            coefficients estimated by polbunch, get the same
                            missing -> misscode treatment applied up front,
                            so both the numest==1 e()-scalar aliases below
                            and the numest>1 report()-scalar matrix columns
                            can just read the "_post" locals directly. This
                            runs unconditionally (not gated on rc==0): when
                            the underlying polbunch call itself failed,
                            elast, se, p, the chi2 stats, and delta_md are
                            still at their missing() initial values from the
                            top of the loop, and the scalar_names bookkeeping
                            below reads these "_post" locals regardless of
                            whether this iteration succeeded -- leaving them
                            undefined on failure previously crashed with
                            "invalid syntax" from a local assignment with
                            nothing on the right-hand side.
                        */
                        local time_post = `time'
                        local p_post    = `p'
                        local se_post   = `se'
                        local rc_post   = `rc'

                        if missing(`time_post') local time_post = `misscode'
                        if missing(`p_post')    local p_post    = `misscode'
                        if missing(`se_post')   local se_post   = `misscode'

                        foreach _tt in minimumdistance hausman {
                            local chi2_`_tt'_post = `chi2_`_tt''
                            local p_`_tt'_post    = `p_`_tt''
                            if missing(`chi2_`_tt'_post') local chi2_`_tt'_post = `misscode'
                            if missing(`p_`_tt'_post')    local p_`_tt'_post    = `misscode'
                        }
                        /*
                            chi2_wald/p_wald are left as genuine Stata
                            missing, not the misscode sentinel, when wald
                            wasn't computed. Unlike hausman/minimum-
                            distance -- which polbunch always attempts
                            for estimator 2/3 -- wald simply isn't part
                            of the default test set for those estimators
                            (only estimator 1/4 run it by default), so
                            missing here means "not applicable" rather
                            than "failed to compute".
                        */
                        local chi2_wald_post = `chi2_wald'
                        local p_wald_post    = `p_wald'
                        local delta_md_post = `delta_md'
                        if missing(`delta_md_post') local delta_md_post = `misscode'

                        if `rc' == 0 {
                            if `numest' == 1 {
                                estimates store `esthold'
                            }
                            else {
                                /*
                                    Real polbunch coefficients (filtered to
                                    report()'s coefficient names when given,
                                    otherwise all of e(b)) and report()-
                                    requested diagnostics both land in the
                                    same combined matrix, under the exact
                                    same [eq]_b[coef] convention: real
                                    coefficients keep polbunch's own equation
                                    names prefixed by the model tag, and
                                    diagnostics ride along as extra
                                    "coefficients" under a <model>_diag
                                    equation. One convention for both means
                                    simulate's default _b collection picks up
                                    everything with no expression list, and
                                    polbunchsim_relabel's reshape needs no
                                    special-casing between the two.
                                */
                                tempname this_b

                                local _want_coefs = ("`report'" == "" | "`_report_coefs'" != "")
                                if `_want_coefs' {
                                    tempname _cof_b
                                    if "`report'" == "" {
                                        matrix `_cof_b' = e(b)
                                    }
                                    else {
                                        mata: _pbsim_filter_b("e(b)", "`_cof_b'", "`_report_coefs'")
                                    }
                                    capture confirm matrix `_cof_b'
                                    if !_rc {
                                        local oldnames : colfullnames `_cof_b'
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
                                                Matrix equation names are
                                                capped at 32 chars. simulate's
                                                own flattened result-variable
                                                name (eq+coef) can still
                                                exceed 32 and fall back to a
                                                generic _sim_N name; don't
                                                fight that here --
                                                polbunchsim_relabel fixes
                                                names up afterward using the
                                                [eq]_b[coef] label simulate
                                                always attaches, regardless
                                                of which name it picked.
                                            */
                                            local neweq = substr("`modelname'_`oldeq'", 1, 32)
                                            local newnames `newnames' `neweq':`oldcoef'
                                        }
                                        matrix colnames `_cof_b' = `newnames'
                                        matrix `this_b' = `_cof_b'

                                        /*
                                            Parallel variance block for a
                                            block-diagonal e(V): this model's
                                            own polbunch covariance on the
                                            diagonal, zeros between models.
                                            Off-model covariance is genuinely
                                            nonzero (same dataset, and btype
                                            only changes the vce, not the point
                                            estimates) -- zeroing it is a
                                            deliberate simplification so that
                                            _se[] is available per model.
                                            Skipped for the whole run if any
                                            contributing model lacks e(V)
                                            (e.g. btype 0 / vce(none)), since a
                                            partial block-diagonal V would be
                                            non-conformable with e(b).
                                        */
                                        capture confirm matrix e(V)
                                        if !_rc {
                                            tempname _cof_V
                                            if "`report'" == "" {
                                                matrix `_cof_V' = e(V)
                                            }
                                            else {
                                                mata: _pbsim_filter_V("e(V)", "`_cof_V'", "`_report_coefs'")
                                            }
                                            capture confirm matrix `_cof_V'
                                            if !_rc & rowsof(`_cof_V') == colsof(`_cof_b') {
                                                matrix rownames `_cof_V' = `newnames'
                                                matrix colnames `_cof_V' = `newnames'
                                                capture confirm matrix `sim_V'
                                                if _rc {
                                                    matrix `sim_V' = `_cof_V'
                                                }
                                                else {
                                                    local _nV1 = rowsof(`sim_V')
                                                    local _nV2 = rowsof(`_cof_V')
                                                    matrix `sim_V' = ///
                                                        ( `sim_V', J(`_nV1', `_nV2', 0) \ ///
                                                          J(`_nV2', `_nV1', 0), `_cof_V' )
                                                }
                                            }
                                            else local _simV_skip = 1
                                        }
                                        else local _simV_skip = 1
                                    }
                                }

                                if "`_report_scals'" != "" {
                                    local _diageq = substr("`modelname'_diag", 1, 32)
                                    foreach _sc of local _report_scals {
                                        local _scval = ``_sc'_post'
                                        /*
                                            chi2_wald/p_wald are left as
                                            genuine missing above (not
                                            applicable, not failed) for the
                                            e()-scalar aliases, but a matrix
                                            posted via -ereturn post- can't
                                            contain missing values ("matrix
                                            has missing values") -- substitute
                                            misscode here, at the point of
                                            building the matrix column, same
                                            as every other diagnostic.
                                        */
                                        if missing(`_scval') local _scval = `misscode'
                                        tempname _scol
                                        matrix `_scol' = (`_scval')
                                        matrix colnames `_scol' = `_sc'
                                        matrix coleq   `_scol' = `_diageq'
                                        capture confirm matrix `this_b'
                                        if _rc matrix `this_b' = `_scol'
                                        else  matrix `this_b' = `this_b', `_scol'
                                    }
                                }

                                capture confirm matrix `this_b'
                                if !_rc matrix `sim_b' = nullmat(`sim_b'), `this_b'
                            }
                        }

                        /*
                            Quantities computed by polbunchsim, rather than
                            coefficients estimated by polbunch, are returned as
                            model-specific e() scalars (the "_post" values were
                            already computed above). simulate can collect them
                            explicitly as e(<model>_<name>).
                        */
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

                e(V) is posted as a block-diagonal matrix -- each model's
                own polbunch covariance on the diagonal, zeros between
                models -- but only when e(b) is purely coefficient columns
                (no report() scalar columns) and every contributing model
                supplied a conformable e(V). Otherwise e(b) is posted
                alone, as before, and per-model elasticity SEs remain
                available in the e(<model>_se) scalars.
            */
            local _post_V = 0
            capture confirm matrix `sim_V'
            if !_rc & "`_simV_skip'" == "" & "`_report_scals'" == "" {
                if rowsof(`sim_V') == colsof(`sim_b') local _post_V = 1
            }
            if `_post_V' {
                * -ereturn post- rejects missing values in V (a failed
                * per-model SE); zero them, matching the off-model blocks.
                mata: st_matrix("`sim_V'", editmissing(st_matrix("`sim_V'"), 0))
                local _bn : colfullnames `sim_b'
                matrix rownames `sim_V' = `_bn'
                matrix colnames `sim_V' = `_bn'
                ereturn post `sim_b' `sim_V'
            }
            else ereturn post `sim_b'

            foreach s of local scalar_names {
                ereturn scalar `s' = ``s''
            }

            ereturn scalar failed = `anyfail'
            ereturn scalar misscode = `misscode'
            ereturn scalar numest = `numest'
            ereturn scalar obs    = `obs'
            ereturn scalar cutoff = `cutoff'
            if `elnum' < .  ereturn scalar el = `elnum'
            else            ereturn local  el "`elasticity'"
            if "`eltrue'" != "" ereturn scalar eltrue = `eltrue'
            if "`ietrue'" != "" ereturn scalar ietrue = `ietrue'
            ereturn local incomeeffect "`incomeeffect'"
            ereturn local buncherror   "`buncherror'"
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
            if "`scalarsonly'" != "" {
                /*
                    Fixed-structure, e(b)-free return for -simulate- /
                    -parallel sim-.  Uses the raw last-iteration locals
                    (genuine missing, not the misscode sentinel).  Every
                    replication returns exactly these scalars, whether or
                    not the underlying polbunch fit converged.
                */
                ereturn clear
                ereturn scalar failed                   = (`final_rc' != 0)
                ereturn scalar sim_rc                    = `final_rc'
                ereturn scalar sim_time                  = `final_time'
                ereturn scalar sim_p                     = `p'
                ereturn scalar sim_se                    = `se'
                ereturn scalar sim_el                    = `elast'
                ereturn scalar sim_b                     = `b_bunch'
                ereturn scalar sim_se_b                  = `se_bunch'
                ereturn scalar sim_p_b                   = `p_bunch'
                ereturn scalar sim_btrue                 = `btrue_use'
                ereturn scalar sim_cover_el             = cond(missing(`p'), ., `p' >= 0.05)
                ereturn scalar sim_cover_b              = cond(missing(`p_bunch'), ., `p_bunch' >= 0.05)
                ereturn scalar sim_chi2_wald             = `chi2_wald'
                ereturn scalar sim_p_wald                = `p_wald'
                ereturn scalar sim_chi2_minimumdistance  = `chi2_minimumdistance'
                ereturn scalar sim_p_minimumdistance     = `p_minimumdistance'
                ereturn scalar sim_chi2_hausman          = `chi2_hausman'
                ereturn scalar sim_p_hausman             = `p_hausman'
                ereturn scalar sim_delta_md              = `delta_md'
                ereturn scalar sim_p_gof                 = `sim_gof_p'
                ereturn scalar sim_gof_deviance          = `sim_gof_deviance'
                ereturn scalar sim_p_gof_below           = `sim_gof_p_below'
                ereturn scalar sim_polyused              = `sim_polyused'
                ereturn scalar sim_perm_rc               = `permrc'
                ereturn scalar sim_perm_p                = `permp'
                ereturn scalar sim_perm_observed         = `permobserved'
                ereturn scalar sim_perm_observed_se      = `permobservedse'
                ereturn scalar sim_perm_n_placebo        = `permnplacebo'
                ereturn scalar sim_perm_n_tried          = `permntried'
                ereturn scalar eltrue                    = `eltrue'
                ereturn scalar ietrue                    = `ietrue'
                ereturn scalar obs                       = `obs'
                ereturn scalar cutoff                    = `cutoff'
                ereturn scalar t0                        = `t0'
                ereturn scalar t1                        = `t1'
                ereturn scalar bw                        = `bw'
                ereturn scalar polynomial                = `polynomial'
                ereturn scalar estimator                 = `estimator'
                ereturn local  cmd "polbunchsim"
                exit
            }
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
                if `elnum' < .  ereturn scalar el = `elnum'
                else            ereturn local  el "`elasticity'"
                if "`eltrue'" != "" ereturn scalar eltrue = `eltrue'
                if "`ietrue'" != "" ereturn scalar ietrue = `ietrue'
                ereturn local incomeeffect "`incomeeffect'"
                ereturn local buncherror   "`buncherror'"
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

                // Scalar columns, under the same <model>_diag convention
                // used for numest>1 (here with no model tag to prefix, since
                // there's only one model): eq = diag, coef = stat name.
                foreach _sc of local _report_scals {
                    local _scval = ``_sc'_post'
                    // See the numest>1 branch: -ereturn post- rejects a
                    // missing value in the matrix even though the
                    // e()-scalar alias is allowed to show genuine missing.
                    if missing(`_scval') local _scval = `misscode'
                    tempname _scol
                    matrix `_scol' = (`_scval')
                    matrix colnames `_scol' = `_sc'
                    matrix coleq   `_scol' = diag
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
            ereturn scalar sim_el              = `elast'
            ereturn scalar sim_b               = `b_bunch'
            ereturn scalar sim_se_b            = `se_bunch'
            ereturn scalar sim_p_b             = `p_bunch'
            ereturn scalar sim_btrue           = `btrue_use'
            ereturn scalar sim_chi2_wald       = `chi2_wald_post'
            ereturn scalar sim_p_wald          = `p_wald_post'
            ereturn scalar sim_chi2_minimumdistance = `chi2_minimumdistance_post'
            ereturn scalar sim_p_minimumdistance    = `p_minimumdistance_post'
            ereturn scalar sim_chi2_hausman    = `chi2_hausman_post'
            ereturn scalar sim_p_hausman       = `p_hausman_post'
            ereturn scalar sim_rc              = `rc_post'

            ereturn scalar obs    = `obs'
            ereturn scalar cutoff = `cutoff'
            if `elnum' < .  ereturn scalar el = `elnum'
            else            ereturn local  el "`elasticity'"
            if "`eltrue'" != "" ereturn scalar eltrue = `eltrue'
            if "`ietrue'" != "" ereturn scalar ietrue = `ietrue'
            ereturn local incomeeffect "`incomeeffect'"
            ereturn local buncherror   "`buncherror'"
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
end