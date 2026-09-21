capture program drop polbunchsim
program polbunchsim, eclass
    syntax [, log reps(integer 1) ///
        obs(integer 5000) cutoff(real 1) ELasticity(string) ///
        t0(real 0.2) t1(real 0.6) bw(real 0.01) ///
        INCOMEeffect(string) buncherror(string) ///
        bootreps(integer 500) POLynomial(numlist integer >=0) ///
        distribution(string) opts(string) ///
        estimator(numlist integer) btype(numlist integer) ///
        clist(string) sample(string) sample4(string)  ///
        est4limits(numlist) limits(numlist) report(string) SCALARsonly DEBUG ///
        PANel(numlist min=3 max=4) B0(real -999) E0(real -999) ///
        HEAP(numlist min=2 max=3) INATTention(real 0) WELfriction(real 0) ///
        PERMStep(real 0) PERMMaxcutoffs(integer 60) NSAMple(integer 0) ///
        FOCAL(numlist min=1 max=3) EXTensive(real 0) TESTBoot(integer 0) TESTBOOTK(numlist integer >=0)]

    /*
        testboot(#) -- default 0 (off).  With # > 1, every estimator 1-4 fit
        is run with polbunch's test(omnibus, boot reps(#)), i.e. an extra
        null-imposed bootstrap of the omnibus test statistic (see polbunch.ado,
        pbx_mdt_bootstat / pbx_wald_bootstat), alongside the usual asymptotic
        chi2 test.  The bootstrap type is polbunch's default, PARAMETRIC
        (fresh multinomial draws around the model's own fitted mean; see
        _pb_parse_testsub); wild/residual are opt-in there and are not
        reachable through testboot().  Its p-value comes back as
        e(sim_p_omnibus_boot) (single fit) or e(<model>_p_omnibus_boot)
        (multi-fit), next to the asymptotic e(sim_p_omnibus)/e(<model>_p_omnibus).
        testbootk(numlist) restricts the bootstrap to those polynomial degrees
        of a polynomial() sweep (implies testboot(199) if testboot() is not
        given); the other degrees still get
        the asymptotic test, from the SAME draw -- so one call gives
        bootstrap and asymptotic p-values on identical data, with no need
        to align random streams across two calls.
        Estimator 0 has no omnibus test and is skipped.  Do not also pass
        test() in opts().  Costs # extra refits per replication.
    */

    /*
        focal(share [width [both]]) / extensive(eps_p) -- passed straight to
        polbunchgendata (see its docstring): excess mass at the cutoff
        without an upstream response, and a participation response above
        the kink.  extensive() drops non-participants, so the raw draw shrinks
        (nsample() is unaffected: it counts rows after sample()).  Neither
        is allowed with panel() where polbunchgendata forbids it.
    */

    /*
        nsample(#) -- default 0 (off).  Requires sample(lo,hi).  Makes the
        number of rows left in the sample() window EXACTLY #, not merely
        obs()*P(window): draw obs() individuals, drop rows outside
        sample(), then keep a uniformly random subset of exactly # of the
        rest (-sample #, count-).  Individuals are iid, so given at least
        # in the window this is an iid draw of # from the observed density
        truncated to the window -- the same distribution as before, with
        the window count fixed rather than binomial.  obs() is then only
        the size of the raw draw: if fewer than # rows land in the window
        the draw is repeated with a larger obs() (see below), and
        e(obs) reports the draw size finally used.  Not with panel().  (Retries:
        obs() grows by 1.25-4x per redraw, at most 8 draws / obs() < 5e7,
        else the replication fails with genrc 459.)
        Without b0(), the realised-bunchers fallback for the coverage
        target is scaled by #/(rows in window before subsampling).
    */

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
        sample(lo,hi) drops rows of the generated data outright (before
        any model is fit), so it applies to every estimator/btype/clist
        combination in the call and shrinks what bootstrap resamples
        from. sample4(lo,hi) mirrors it but only as an `if inrange(z,
        lo,hi)' qualifier on the estimator-4 (Saez) polbunch call --
        the underlying dataset (and every other estimator) is
        untouched. Use sample() to change what all estimators see;
        use sample4() on top of it to additionally narrow (or widen)
        just the Saez window, since Saez's linear-restriction estimator
        is typically fit on a narrower symmetric window than the
        polynomial estimators.
    */

    /*
        polynomial() takes a numlist, same as btype()/estimator(): every
        requested degree is fit SEPARATELY on the one draw generated for
        this replication (like sweeping btype()/estimator()/clist(), not
        like sample()'s once-only restriction), and the model tag gets a
        "k<degree>" component (e.g. "e3k2") whenever more than one degree
        is requested -- a single value adds no tag, same as every other
        dimension. Estimator 4 (Saez) ignores its own polynomial() value
        internally (polbunch always fits a two-point flat counterfactual
        there), so sweeping k for estimator 4 alone is a harmless no-op.
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
        if "`btype'" == "" local btype 1
        if "`estimator'" == "" local estimator 3
        if "`polynomial'" == "" local polynomial 1
        if "`elasticity'" == "" local elasticity 0.4
        local elnum = real("`elasticity'")

        if "`clist'" == "" local clist `"noconstant"'

        local misscode = -1e300

        if `nsample' < 0 {
            noi di as error "nsample() must be a positive integer"
            exit 198
        }
        if `nsample' > 0 & "`sample'" == "" {
            noi di as error "nsample() requires sample(lo,hi)"
            exit 198
        }
        if "`testbootk'" != "" & `testboot' <= 1 local testboot 199
        if `nsample' > 0 & "`panel'" != "" {
            noi di as error "nsample() is not allowed with panel()"
            exit 198
        }

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
        local _scalarnames " chi2_wald p_wald chi2_hausman p_hausman chi2_minimumdistance p_minimumdistance chi2_omnibus p_omnibus delta_md p_gof deviance p_gof_below polyused time se p rc "
        if "`report'" != "" {
            foreach _item of local report {
                if inlist("`_item'", "hausman", "wald", "minimumdistance", "omnibus") {
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
        local numk : word count `polynomial'
        local numest = `numb' * `nume' * `numc' * `numk'

        /*
            Only stamp a dimension into the model tag when it actually
            varies in this call, so e.g. a run that only sweeps
            estimator()/clist() gets tags like "e3c0" instead of the
            noisier "b1e3c0" -- and a run that only sweeps btype() gets
            "b1", "b2", etc. Same for polynomial(): a numlist sweep tags
            "k2", "k3", ...; a single value (the common case) adds no tag.
        */
        local vary_bt = (`numb' > 1)
        local vary_e  = (`nume' > 1)
        local vary_c  = (`numc' > 1)
        local vary_k  = (`numk' > 1)

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
        if "`focal'" != "" local heapopt `heapopt' focal(`focal')
        if `extensive' > 0 local heapopt `heapopt' extensive(`extensive')

        // draw; with nsample(), redraw larger until >= nsample rows fall in sample()
        local _try  = 0
        local _done = 0
        local _nwin = .
        while !`_done' {
            local ++_try
            local _done = 1
            clear
            capture noisily polbunchgendata z, obs(`obs') cutoff(`cutoff') ///
                elasticity(`elasticity') t0(`t0') t1(`t1') `log' distribution(`distribution') ///
                incomeeffect(`incomeeffect') buncherror(`buncherror') `panelopt' `heapopt' ///
                inattention(`inattention') welfriction(`welfriction')

            local genrc = _rc
            if `genrc' == 0 {
                local _r_el = r(el_mean)
                local _r_ie = r(incomeeffect)
                local _r_nb = r(n_bunchers)
                if `nsample' > 0 {
                    count if inrange(z, `sample')
                    local _nwin = r(N)
                    if `_nwin' < `nsample' {
                        if `_try' < 8 & `obs' < 5e7 {
                            local obs = ceil(`obs' * min(4, max(1.25, 1.1*`nsample'/max(`_nwin',1))))
                            local _done = 0
                        }
                        else {
                            noi di as error "nsample(`nsample'): only `_nwin' rows in sample() after `_try' draws (obs=`obs')"
                            local genrc = 459
                        }
                    }
                }
            }
        }

        local btrue_realized = .
        if `genrc' == 0 {
            local eltrue = cond(`e0' >= 0, `e0', `_r_el')
            local ietrue = `_r_ie'
            local btrue_realized = `_r_nb'
            if `nsample' > 0 local btrue_realized = `_r_nb' * `nsample' / `_nwin'
        }
        if "`eltrue'" == "" local eltrue = cond(`e0' >= 0, `e0', `elnum')
        if "`eltrue'" == "" local eltrue = .
        if "`ietrue'" == "" local ietrue = .

        // population excess-mass count for the coverage test
        local btrue_use = cond(`b0' >= 0, `b0', `btrue_realized')

        if `genrc' == 0 {
            if "`sample'" != "" {
                drop if !inrange(z, `sample')
                if `nsample' > 0 sample `nsample', count
            }

            foreach bt of numlist `btype' {
                foreach e of numlist `estimator' {
                    foreach c in `clist' {
                        foreach k of numlist `polynomial' {

                        if "`c'" == "constant" local cval = 1
                        else local cval = 0

                        local modelname ""
                        if `vary_bt' local modelname "`modelname'b`bt'"
                        if `vary_e'  local modelname "`modelname'e`e'"
                        if `vary_c'  local modelname "`modelname'c`cval'"
                        if `vary_k'  local modelname "`modelname'k`k'"

                        timer clear
                        timer on 1

                        if `e' == 4 & "`sample4'" != "" local iff `"if inrange(z, `sample4')"'
                        else local iff

                        if `e' == 4 & "`est4limits'" != "" {
                            local uselimits limits(`est4limits')
                        }
                        else {
                            local uselimits limits(`limits')
                        }

                        local optsx `opts'
                        if `testboot' > 1 & `e' != 0 {
                            local _doboot 1
                            if "`testbootk'" != "" local _doboot : list k in testbootk
                            if `_doboot' local optsx `opts' test(omnibus, boot reps(`testboot'))
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
                                pol(`k') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') vce(none) `c' ///
                                 `uselimits' `optsx'
                        }
                        else if `bt' == 1 {
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`k') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') vce(analytic) `c' ///
                                `optsx' `uselimits'
                        }
                        else if `bt' == 2 {
                            /* individual (person-year) nonparametric bootstrap,
                               Saez-2010 style: resample rows, re-bin, re-fit */
                            capture `dbgpfx' bootstrap ///
                                _bs_el = _b[bunching:elasticity] ///
                                _bs_nb = _b[bunching:number_bunchers], ///
                                reps(`bootreps') nodots: ///
                                polbunch z `iff', cutoff(`cutoff') ///
                                pol(`k') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') vce(none) `c' ///
                                `optsx' `uselimits'
                        }
                        else if `bt' == 3 {
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`k') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') vce(bootstrap) ///
                                 `c' `optsx' `uselimits'
                        }
                        else if `bt' == 4 {
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`k') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') vce(bayes) ///
                                 `c' `optsx' `uselimits'
                        }
                        else if `bt' == 5 {
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`k') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') vce(bootstrap) ///
                                `c' `optsx' `uselimits' nozero
                        }
                        else if `bt' == 6 {
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`k') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') vce(bayes) ///
                                nozero `c' `optsx' `uselimits'
                        }
                        else if `bt' == 7 {
                            /* analytic Eicker-White sandwich */
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`k') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') vce(robust) `c' ///
                                `optsx' `uselimits'
                        }
                        else if `bt' == 8 {
                            /* analytic cluster-robust (needs panel -> pid) */
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`k') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') vce(cluster pid) `c' ///
                                `optsx' `uselimits'
                        }
                        else if `bt' == 9 {
                            /* Chetty/CFOP residual bootstrap over bins */
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`k') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') bootreps(`bootreps') ///
                                vce(bootstrap, residual) `c' `optsx' `uselimits'
                        }
                        else if `bt' == 10 {
                            /* person-clustered nonparametric bootstrap (benchmark
                               for vce(cluster) under pooled-panel dependence) */
                            capture `dbgpfx' bootstrap ///
                                _bs_el = _b[bunching:elasticity] ///
                                _bs_nb = _b[bunching:number_bunchers], ///
                                reps(`bootreps') cluster(pid) nodots: ///
                                polbunch z `iff', cutoff(`cutoff') ///
                                pol(`k') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') vce(none) `c' ///
                                `optsx' `uselimits'
                        }
                        else if `bt' == 11 {
                            /* quasi-multinomial: conventional V scaled by the
                               Pearson dispersion phi-hat (Dirichlet-multinomial
                               / Kish design-effect fallback when M is
                               unavailable -- histogram-only) */
                            capture `dbgpfx' polbunch z `iff', cutoff(`cutoff') ///
                                pol(`k') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') vce(conventional) scale(x2) `c' ///
                                `optsx' `uselimits'
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
                                pol(`k') bw(`bw') t0(`t0') t1(`t1') ///
                                `log' estimator(`e') vce(analytic) `c' ///
                                `optsx' `uselimits'
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
                        local p_omnibus_boot = .
                        local omnibus_boot_n = .
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

                            /*
                                Omnibus test: polbunch itself now reports
                                its one default nested-restriction test per
                                estimator (wald for 1/4, minimumdistance for
                                2/3) under e(chi2_omnibus)/e(p_omnibus), so
                                just read it straight off -- no need to
                                pick between chi2_wald/chi2_minimumdistance
                                by estimator number here any more. Missing
                                for estimator 0 (no nested restriction to
                                test) and for an off-label diagnostic test
                                requested via opts(test(wald)) on estimator
                                2/3 (which polbunch deliberately keeps under
                                its own chi2_wald/p_wald, not omnibus, since
                                it is NOT that estimator's own default test
                                -- see polbunch.ado's test() docs).
                            */
                            capture local chi2_omnibus = e(chi2_omnibus)
                            if _rc local chi2_omnibus = .
                            capture local p_omnibus = e(p_omnibus)
                            if _rc local p_omnibus = .
                            capture local p_omnibus_boot = e(p_omnibus_boot)
                            if _rc local p_omnibus_boot = .
                            capture local omnibus_boot_n = e(omnibus_boot_n)
                            if _rc local omnibus_boot_n = .

                            /*
                                Backfill the legacy chi2_wald/p_wald and
                                chi2_minimumdistance/p_minimumdistance
                                locals from the new e(chi2_omnibus) for the
                                NORMAL (non-off-label) case, so existing
                                callers reading polbunchsim's own
                                sim_chi2_wald/sim_p_wald (estimator 1/4) or
                                sim_chi2_minimumdistance/sim_p_minimumdistance
                                (estimator 2/3) -- e.g. sim_spectest.do,
                                spworker.do -- keep working unchanged. Only
                                fills in when the raw local came back
                                missing (i.e. NOT the off-label case, where
                                polbunch itself still populates chi2_wald
                                directly and that real value must win).
                            */
                            if inlist(`e', 1, 4) & missing(`chi2_wald') {
                                local chi2_wald = `chi2_omnibus'
                                local p_wald    = `p_omnibus'
                            }
                            if inlist(`e', 2, 3) & missing(`chi2_minimumdistance') {
                                local chi2_minimumdistance = `chi2_omnibus'
                                local p_minimumdistance    = `p_omnibus'
                            }

                            /*
                                Reference-region GoF/deviance and the
                                actually-fitted polynomial order (polbunch
                                silently lowers it by 1+ on rank
                                deficiency -- see its own "Note: Polynomial
                                order lowered..."). These are set by
                                polbunch's MAIN fit itself (ereturn scalar
                                deviance_p / polynomial / ...) regardless
                                of vce()/btype, so no `estat gof' call is
                                needed -- read them directly off e().
                                Computed for EVERY estimator, not just 0:
                                for estimator 0 (no nested restriction to
                                run an omnibus test against) this is the
                                ONLY assumption test available; for
                                estimators 1-3 it's the pooled-deviance
                                view of the same overall null the omnibus
                                test targets, via the many-bin reference
                                fit rather than the structural contrast.
                                Not populated after a bootstrap-PREFIX
                                btype (2, 10): the wrapping `bootstrap'
                                command's own ereturn post does not retain
                                the wrapped command's extra scalars.
                            */
                            capture local sim_polyused = e(polynomial)
                            if _rc local sim_polyused = .
                            capture local sim_gof_p = e(deviance_p)
                            if _rc local sim_gof_p = .
                            capture local sim_gof_deviance = e(deviance)
                            if _rc local sim_gof_deviance = .
                            capture local sim_gof_p_below = e(deviance_below_p)
                            if _rc local sim_gof_p_below = .
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

                        /*
                            chi2_omnibus/p_omnibus follow the same
                            "genuine missing == not applicable" convention
                            as chi2_wald/p_wald above (estimator 0 has no
                            omnibus test). GoF/polyused are genuinely
                            missing only on failure to compute, never
                            "not applicable" (GoF is defined for every
                            estimator) -- same convention either way,
                            since the report()-scalar matrix-column loop
                            below substitutes misscode at the point of
                            use regardless.
                        */
                        local chi2_omnibus_post = `chi2_omnibus'
                        local p_omnibus_post    = `p_omnibus'
                        local p_omnibus_boot_post = `p_omnibus_boot'
                        local omnibus_boot_n_post = `omnibus_boot_n'
                        local p_gof_post        = `sim_gof_p'
                        local deviance_post     = `sim_gof_deviance'
                        local p_gof_below_post  = `sim_gof_p_below'
                        local polyused_post     = `sim_polyused'

                        if `rc' == 0 {
                            if `numest' == 1 {
                                estimates store `esthold'
                            }
                            else if "`scalarsonly'" == "" {
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
                            `modelname'_chi2_omnibus `modelname'_p_omnibus ///
                            `modelname'_p_omnibus_boot ///
                            `modelname'_omnibus_boot_n ///
                            `modelname'_chi2_hausman `modelname'_p_hausman ///
                            `modelname'_p_gof `modelname'_deviance ///
                            `modelname'_p_gof_below `modelname'_polyused ///
                            `modelname'_rc

                        local `modelname'_time                 = `time_post'
                        local `modelname'_se                   = `se_post'
                        local `modelname'_p                    = `p_post'
                        local `modelname'_chi2_omnibus         = `chi2_omnibus_post'
                        local `modelname'_p_omnibus            = `p_omnibus_post'
                        local `modelname'_p_omnibus_boot       = `p_omnibus_boot_post'
                        local `modelname'_omnibus_boot_n       = `omnibus_boot_n_post'
                        local `modelname'_chi2_hausman         = `chi2_hausman_post'
                        local `modelname'_p_hausman            = `p_hausman_post'
                        local `modelname'_p_gof                = `p_gof_post'
                        local `modelname'_deviance             = `deviance_post'
                        local `modelname'_p_gof_below          = `p_gof_below_post'
                        local `modelname'_polyused             = `polyused_post'
                        local `modelname'_rc                   = `rc_post'

                        if `rc' local anyfail = 1

                        local final_bt   `bt'
                        local final_e    `e'
                        local final_cval `cval'
                        local final_k    `k'
                        local final_rc   `rc'
                        local final_time `time'
                        }
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
            if "`scalarsonly'" != "" {
                /*
                    Fixed-structure, e(b)-free return, same spirit as the
                    numest==1 scalarsonly branch below: no coefficient
                    matrix is built at all (so one model's failure never
                    shrinks e(b)'s column count and breaks -simulate-'s
                    per-replication structure), just the per-model
                    e(<model>_*) scalars set above, always present with
                    the same names whether or not each model's fit
                    succeeded.
                */
                ereturn clear
                foreach s of local scalar_names {
                    ereturn scalar `s' = ``s''
                }
                ereturn scalar failed  = `anyfail'
                ereturn scalar misscode = `misscode'
                ereturn scalar numest  = `numest'
                ereturn scalar obs     = `obs'
                ereturn scalar cutoff  = `cutoff'
                if `elnum' < .  ereturn scalar el = `elnum'
                else            ereturn local  el "`elasticity'"
                if "`eltrue'" != "" ereturn scalar eltrue = `eltrue'
                if "`ietrue'" != "" ereturn scalar ietrue = `ietrue'
                ereturn local incomeeffect "`incomeeffect'"
                ereturn local buncherror   "`buncherror'"
                ereturn scalar t0     = `t0'
                ereturn scalar t1     = `t1'
                ereturn scalar bw     = `bw'
                ereturn local polynomial "`polynomial'"
                ereturn scalar bootreps   = `bootreps'
                ereturn local btype "`btype'"
                ereturn local estimator "`estimator'"
                ereturn local clist "`clist'"
                ereturn local cmd "polbunchsim"
                exit
            }

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
            ereturn local polynomial "`polynomial'"
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
                ereturn scalar sim_chi2_omnibus          = `chi2_omnibus'
                ereturn scalar sim_p_omnibus             = `p_omnibus'
                ereturn scalar sim_p_omnibus_boot        = `p_omnibus_boot'
                ereturn scalar sim_omnibus_boot_n        = `omnibus_boot_n'
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
                ereturn local polynomial "`polynomial'"
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
            ereturn scalar sim_chi2_omnibus    = `chi2_omnibus_post'
            ereturn scalar sim_p_omnibus       = `p_omnibus_post'
            ereturn scalar sim_p_omnibus_boot  = `p_omnibus_boot_post'
            ereturn scalar sim_omnibus_boot_n   = `omnibus_boot_n_post'
            ereturn scalar sim_chi2_hausman    = `chi2_hausman_post'
            ereturn scalar sim_p_hausman       = `p_hausman_post'
            ereturn scalar sim_p_gof           = `p_gof_post'
            ereturn scalar sim_gof_deviance    = `deviance_post'
            ereturn scalar sim_p_gof_below     = `p_gof_below_post'
            ereturn scalar sim_polyused        = `polyused_post'
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