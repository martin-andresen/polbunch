*! polbunchbias -- analytical bias of polynomial bunching estimators
*! version 2.5.0  07sep2026
*!
*! 2.5.0: estimator(4) + constant is now refused.  Saez with a
*!        constant-density inversion applies a FLAT counterfactual in both
*!        the excess-mass and the response-inversion step, so the two
*!        stages agree and the internal-consistency bias this command
*!        reports is zero by construction (the same situation as
*!        estimators 0/3 under exact+splitmass).  This is not a statement
*!        that Saez/constant is unbiased -- it is biased whenever the
*!        counterfactual slopes over the excluded window -- so the command
*!        errors rather than print a deterministic zero that reads as
*!        "unbiased".  estimator(4) + exact keeps its non-trivial internal
*!        check (step-function mass vs. the implied linear counterfactual).
*!        polbunch's inline-bias call already swallows this error silently,
*!        exactly as it does for estimators 0/3.  Also: both the standalone
*!        and the inline bias table now carry a note that the figure is
*!        measured against that estimator's own counterfactual and is not
*!        comparable across estimators.
*!
*! 2.4.0: estimator-2 + splitmass bias_B HARMONISED with polbunch.  polbunch
*!        forms the splitmass Bhat as Hstar - int_{zL}^{z*} h0hat -
*!        int_{z*}^{zH} h0hat/(1+deltahat) -- it deltahat-deflates the
*!        FITTED h0hat above the kink (the Chetty h1 = h0/(1+delta)
*!        restriction).  This file instead subtracted an iso-elastic
*!        (1 - 1/x) * int_0^H h0_TRUE term, which predicted a ~-0.03
*!        elasticity bias polbunch does not actually exhibit (an oracle
*!        plim check with the true h0 put the splitmass fixed point 0.03
*!        off the truth).  The extra term over poolmass is now
*!        (delc/(1+delc)) * int_{z*}^{zH} h0hat / bw with delc the NLS
*!        pseudo-true delta and h0hat = betat + fitting-stage bias -- an
*!        exact algebraic consequence of following polbunch's B formula at
*!        any degree.  Oracle splitmass gap 0.03 -> 0.002.  polbunchbias_legacy.do
*!        (K=1 test reference) and test_pbx_est2.do B1 updated to match.
*!
*! 2.4.0: `iterate' now solves the estimator-1 and estimator-2
*!        self-consistency fixed point with a single-start Newton
*!        (finite-difference Jacobian + backtracking line search;
*!        pbx_bias_newton_ws / pbx_bias_scres_ws), started from the fitted
*!        h0 shape and the naive elasticity.  For estimator 2 the damped
*!        fixed-point substitution is not a contraction and ran away in
*!        ~100% of degree->=5 draws (forcing the un-iterated /
*!        lower-order-promotion fallbacks to carry the correction); for
*!        estimator 1 it converges only slowly.  The fixed point is unique
*!        and sits at the truth -- oracle plim check, both mass axes for
*!        est 2 (needed the splitmass bias_B harmonisation above first) --
*!        and Newton reaches it directly: MC RMSE ~0.02, nearly unbiased,
*!        0 fallbacks across degree-3/5 counterfactuals, level and log,
*!        a range of elasticities/windows.  The damped loop is unchanged
*!        and stays as the fallback when Newton does not converge to an
*!        admissible root, and remains the only path for estimators 0/3/4
*!        (estimator 3's damped `iterate' already reaches the truth -- it
*!        has no fitting-stage circularity).  Bit-identical to 2.3.0
*!        whenever `iterate' is not requested or the estimator is 0/3/4.
*!        KNOWN RESIDUAL: estimator 2 + log mode on a low-order
*!        counterfactual leaves a ~0.04 elasticity bias (pre-existing in
*!        the poolmass path, still far better than the old ~naive
*!        fallback); estimator 1 + log and estimator 3 + log are clean.
*!
*! 2.3.0: also -- the estimator-2 (K+2)-square solve's diagonal
*!        preconditioner now scales its delta index by 1/sqrt(qff) instead
*!        of 1.  The delta row/column scales with the overall HEIGHT of h0
*!        (Mhi linear in the coefficients, qff = b*Ghi*b' quadratic) while
*!        the beta block does not, so a tall h0 -- counts ~1e5, or an
*!        h0poly() supplied in count units -- used to make GtG[np,np]
*!        dominate and trip the conditioning guard ("counterfactual not
*!        identified ... at any order down to 1") even on a perfectly
*!        solvable system.  Exact change of basis; bit-identical to before
*!        on normal e()-mode input (whose small normalised window already
*!        keeps qff in range), only unblocks large-magnitude standalone
*!        h0poly() calls.
*!
*! 2.3.0: estimator 2's delta-solve is now multi-start (pbx_bias_e2_solve /
*!        pbx_bias_e2_fixedpoint), replacing the single Gauss-Newton path
*!        started at delc=Delta that 2.2.0 flagged as needing this fix --
*!        the window-width sweep that motivated it showed wild, sign-
*!        flipping bias_e swings across neighbouring window widths, all
*!        traced to that single trajectory wandering through a badly-
*!        conditioned region as delc moves away from its start.  The FIXED
*!        POINT EQUATION itself is unchanged (same (K+2)-square GtG/Gtu
*!        assembly, same Mhi/qff Jacobian terms built from the fixed true
*!        coefficients, same 0.7-damped step, same boundary clamp, same
*!        post-loop conditioning re-check) -- only the starting value is
*!        now a grid of ~10 points instead of one, with ties among multiple
*!        converged roots broken by a concentrated-SSR merit function
*!        (pbx_bias_e2_objQ).  An earlier draft of this fix instead
*!        reformulated the equation itself (profiling beta out and globally
*!        minimizing that same concentrated SSR via golden section, with no
*!        multi-start at all) -- mathematically a well-posed problem, but
*!        validating it against test_pbx_est2.do's independent K=1
*!        reference (pb_fobias_core) showed differences up to 0.10: it
*!        silently changed the estimand, not just its robustness, because
*!        the Gauss-Newton Jacobian being replaced holds the true
*!        coefficients fixed rather than re-differentiating at the running
*!        bias estimate, which is a deliberate (already-validated) choice,
*!        not an approximation error.  The multi-start-of-the-original-
*!        equation design here reproduces that reference to ~1e-12 on the
*!        full K=1 grid (previously machine precision, ~1e-15, against the
*!        single-start method it replaces) and 0/8 failures on the K>=2
*!        closed-form-vs-brute-force suite, while genuinely recovering
*!        additional identified window widths in the chaotic region that
*!        the single start missed, and correctly reporting missing (rather
*!        than a plausible-looking wrong number) where no starting value
*!        converges to a well-conditioned root.  See
*!        polbunchbias-exact-rewrite.md (memory) for the full derivation,
*!        including the two dead ends hit before landing on this design.
*!
*! 2.2.0: the bias engine (pbx_bias_core/pbx_bias_solve) now works NATIVELY
*!        in polbunch's own normalised coordinate w = (z - z*)/xscale
*!        instead of converting to raw z-units up front.  Previously,
*!        e()-mode read the fitted h0 coefficients from e(b) (fit on
*!        polbunch's normalised regressor) and immediately shifted AND
*!        rescaled them into raw z-units (pbx_ecoef_transform) before doing
*!        any bias arithmetic; the rescale-by-1/xscale^k step stretches the
*!        coefficient vector across many orders of magnitude whenever
*!        xscale is far from 1 (a wide window or a cutoff far from 0), and
*!        that damage was done before the analytical-bias machinery (Gram
*!        solves, the iterate fixed point, the h0poly/relslope fallbacks)
*!        ever saw the polynomial.  Now the shift (still the numerically
*!        delicate step, still synthetic division) happens ALONE, staying
*!        in w-space; the rescale is never done at all -- z-unit-only
*!        quantities (bias_response, bias_slope, biasbeta, the zstar/window
*!        echo columns) are converted back to z-units in ONE place, at the
*!        very end, instead of up front on the input.  This required
*!        decoupling zstar's old dual role (window-centring POSITION, now
*!        simply 0 since w is already centred at the cutoff, vs. the
*!        Delta-to-level-shift ECONOMIC MAGNITUDE, now zstarw = zstar/
*!        xscale) and the equivalent rho (log mode's analogous shift)
*!        dual role.
*!
*!        Developed and validated as a separate file (polbunchbias_norm.ado,
*!        pbn_* Mata namespace) before being promoted here, specifically so
*!        it could be tested side-by-side against this file with zero risk
*!        of the two colliding.  Validated: all 5 estimators bit-identical
*!        (level mode) to the pre-2.2.0 engine except estimator 2's own NLS
*!        delta-solve, which agrees to ~4e-6 relative (a different but
*!        equally valid floating-point path through an iterative solve, not
*!        a formula difference); log mode validated bit-identical for
*!        estimators 1/3/4 (2/0 not reachable in the validation DGP used --
*!        hasresp=0/weak identification there -- but share the identical
*!        fix pattern already confirmed in 1/3/4).  A genuine log-mode unit
*!        bug (rho and the relative-slope prefactor both needed the same
*!        w-space correction as zstar; missed on the first pass since only
*!        level mode was tested) was caught by this validation and fixed
*!        before promotion -- see polbunchbias-exact-rewrite.md (memory)
*!        for the full derivation.  Measured benefit: modest, not dramatic
*!        (a direct round-trip shift-precision test showed ~20-30% error
*!        reduction at an extreme synthetic shift, most of the damage
*!        having already been mitigated by the iterate-loop preconditioner
*!        added earlier); preferred anyway since it is a strict correctness
*!        improvement with no measured downside.  Does NOT address the
*!        separate estimator-2 NLS chaotic-sensitivity issue (wild swings
*!        in the bias across neighbouring window widths) -- that is
*!        intrinsic to the Gauss-Newton delta-solve, present identically in
*!        both coordinate systems, and needs its own (multi-start-style) fix.
*!
*! 2.1.0: estimators 1/2 -- plausibility guard on the fitting-stage bias
*!        (a coefficient bias exceeding the counterfactual it corrects is a
*!        near-collinear-design artifact, not an economic bias) and an
*!        automatic polynomial-order fallback: on a wide / off-centre window
*!        at high order the bias is recomputed with the best lower-degree L2
*!        approximation of the counterfactual, reporting the order used and
*!        the L2 change (r(bias_polynomial), r(bias_polyfull), r(bias_l2err)).
*!        New nofallback option.  Also: when `iterate' fails to converge at
*!        the reported order, search for the highest lower order where it
*!        DOES reach a genuine self-consistent fixed point.  A hand-fitted
*!        or digitized h0 is noisiest at high order, and the bias integral
*!        extrapolates h0 across the excluded window -- exactly where that
*!        noise is amplified -- so BY DEFAULT the converged lower-order fit
*!        is PROMOTED to the primary reported bias, trading a small
*!        deliberate truncation for materially less extrapolation variance
*!        (r(bias_uniter_*) keeps the demoted full-order plug-in number,
*!        never discarded).  New nopromote option restores the pre-04sep2026
*!        default (full-order plug-in primary, converged fit in
*!        r(bias_iter_*)).  Promotion is gated by promotetol(#) (default
*!        0.3): only promote if the L2 shape change to reach convergence is
*!        below this -- MC validation (known-truth simulation) found that
*!        forcing convergence can require collapsing h0 almost entirely
*!        (e.g. estimator 2 often only converges at K=1), which is a WORSE
*!        bias correction than the un-iterated fit, not a better one; a
*!        modest L2 change (as typical for estimator 1) IS an improvement.
*!
*!        iterate's own arithmetic is now done in a diagonally
*!        preconditioned basis (same change of basis as the Gram solves)
*!        instead of on raw cutoff-centred coefficients directly -- those
*!        can span many orders of magnitude purely from polbunch's own
*!        internal normalisation, a source of spurious divergence distinct
*!        from genuine model misspecification.  New allownegative option
*!        (default off, i.e. floored): the corrected elasticity is floored
*!        at 0 by default (Slutsky: a compensated elasticity is
*!        non-negative), mirroring polbunch's delta>=0 default for the
*!        convex-kink profile -- validated to never bind in testing, so it
*!        costs nothing when it isn't needed.  New h0check option (default
*!        OFF, i.e. NOT checked): optionally also reject a candidate whose
*!        implied counterfactual goes negative anywhere in the fitting
*!        region.  Tried as default-on and found, via MC validation against
*!        a known true bias, to erase the entire benefit of promotion
*!        rather than just filter bad candidates -- h0 here is a low-degree
*!        polynomial APPROXIMATION to a noisily estimated density, and such
*!        an approximation can legitimately dip slightly negative in a
*!        low-density region without the true h0 being invalid; unlike the
*!        elasticity floor this is a threshold on an entire noisy function
*!        shape, not a single theoretically-bounded scalar.
*!

*! The counterfactual density h0 is treated as an arbitrary degree-K
*! polynomial in the centred running variable s (s = z - z* in levels,
*! s = ln z - ln z* in logs).  For a polynomial DGP the reported biases
*! are EXACT, not a first-order approximation; a degree-1 polynomial (a
*! constant and a slope) recovers the earlier local-linear behaviour.
*!
*! Modes:
*!   e()      -- after polbunch: h0 (degree e(polynomial)) is read from
*!               e(b) and recentred on the cutoff.  log mode assumes the
*!               running variable is already ln(earnings), as in polbunch.
*!               The response-inversion axis (constant/exact) and the
*!               bunching-mass axis (poolmass/splitmass) are taken from
*!               e(constant)/e(nosplit) unless overridden on the command.
*!   explicit -- 9 primary options plus relslope() (degree 1) or h0poly()
*!               (the centred, height-first coefficient list, any degree).
*!
*! The model-implied-B option (bmodel) has been removed: estimator 2 always
*! reports the reduced-form bunching mass, governed by poolmass/splitmass.
*!
*! Implementation: pbx_bias_core() is the arbitrary-degree bias engine;
*! pbx_bias_solve() adds the optional self-consistent iteration; the
*! pbx_*() Mata layer below is a small polynomial-algebra toolkit.

capture program drop polbunchbias
program define polbunchbias, rclass
    version 16.0

    syntax [, ESTimator(numlist max=1) ZSTAR(numlist max=1) ///
        T0(numlist max=1) T1(numlist max=1) ///
        RELSLOPE(numlist max=1) H0poly(numlist) ELasticity(numlist max=1) ///
        ZLO(numlist max=1) ZHI(numlist max=1) ///
        ZL(numlist max=1) ZH(numlist max=1) ///
        LOG ///
        BW(numlist max=1) ITERate TOLerance(real 1e-10) ///
        MAXITER(integer 1000) UNDERRELAX(real 0.5) ///
        CONstant EXACT SPLITmass POOLmass NOFALLback NOPROMOTE ///
        PROMOTEtol(real 0.3) ALLOWNEGATIVE H0check ]

    // vestigial: the model-implied-B option was removed (mirrors polbunch's
    // dropped Bmodel).  Estimator 2 always uses the reduced-form bunching
    // mass now, governed by poolmass / splitmass.  Kept as an always-0
    // argument so the bias engine's 26-column return layout is unchanged.
    local bmodel = 0

    if "`relslope'" != "" & "`h0poly'" != "" {
        di as error "polbunchbias: specify at most one of relslope() / h0poly()"
        exit 198
    }
    if "`constant'" != "" & "`exact'" != "" {
        di as error "polbunchbias: specify at most one of constant / exact"
        exit 198
    }
    if "`poolmass'" != "" & "`splitmass'" != "" {
        di as error "polbunchbias: specify at most one of poolmass / splitmass"
        exit 198
    }

    // count of the 9 primary options (relslope/h0poly do NOT participate --
    // "polbunchbias, relslope(#)" must still land in e() mode)
    local ncore = ("`estimator'"!="") + ("`zstar'"!="") + ///
        ("`t0'"!="") + ("`t1'"!="") + ("`elasticity'"!="") + ///
        ("`zlo'"!="") + ("`zhi'"!="") + ("`zl'"!="") + ("`zh'"!="")

    tempname craw bcoef
    local _e_nosplit
    local _e_constant
    // xscale: "what units is bcoef expressed in, relative to raw z?"  Only
    // the auto-extracted e(b) polynomial coefficients are in polbunch's own
    // normalised units (xscale = e(xscale)); h0poly()/relslope()/the Saez
    // two-point construction are always user-supplied or built directly in
    // raw z (or log-earnings) units, same as explicit mode -- xscale=1 for
    // all of those, overridden only in the one branch below that reads
    // e(xscale).  pbx_bias_core/pbx_bias_solve work entirely in this
    // xscale's units; xscale=1 recovers exactly the old raw-coordinate
    // behaviour, which is the main regression check for this rewrite.
    local _xscale = 1

    if `ncore' == 0 {
        // ================= e() MODE =================
        if "`e(cmd)'" != "polbunch" {
            di as error "No options supplied, but last estimation command was not polbunch"
            if "`e(cmd)'" != "" di as error "(found: e(cmd) = `e(cmd)')"
            exit 301
        }
        if "`e(transform)'" == "notransform" {
            di as error "polbunchbias requires transformed polbunch output"
            di as error "(re-run polbunch without the notransform option)"
            exit 321
        }

        local estimator = e(estimator)
        local zstar     = e(cutoff_orig)
        local _bw       = e(bw_orig)
        local zl        = e(lower_limit)
        local zh        = e(upper_limit)
        local islog_val = e(log)
        local zlo       = e(zlo)
        local zhi       = e(zhi)
        local _zname    = e(zname)

        capture confirm scalar e(nosplit)
        if !_rc local _e_nosplit = e(nosplit)
        capture confirm scalar e(constant)
        if !_rc local _e_constant = e(constant)

        if missing(`zlo') | missing(`zhi') {
            di as error "e(zlo)/e(zhi) not found; re-run polbunch (updated version required)"
            exit 111
        }
        local _t0 = e(t0)
        local _t1 = e(t1)
        if missing(`_t0') | missing(`_t1') {
            di as error "e(t0)/e(t1) not found; re-run polbunch with t0() and t1() specified"
            exit 111
        }
        local t0 = `_t0'
        local t1 = `_t1'

        capture local _elast = _b[bunching:elasticity]
        if _rc | missing(`_elast') {
            di as error "bunching:elasticity not found in e(b)"
            di as error "(polbunch must be run with t0() and t1())"
            exit 111
        }
        local elasticity = `_elast'

        if "`bw'" == "" local bw = `_bw'

        // ---- build the counterfactual polynomial bcoef ----
        if "`h0poly'" != "" {
            local _nb : word count `h0poly'
            matrix `bcoef' = J(1, `_nb', 0)
            local _jb = 0
            foreach _v of local h0poly {
                local ++_jb
                matrix `bcoef'[1,`_jb'] = `_v'
            }
        }
        else if "`relslope'" != "" {
            // relslope is the relative slope w.r.t. the running variable
            // (per unit z in levels, per unit ln z in logs)
            matrix `bcoef' = (1, `relslope')
        }
        else if `estimator' == 4 {
            // Saez: recover the slope from the implicit two-point counterfactual
            local _hminus = _b[h0:_cons]
            local _hplus  = _b[h1:_cons]
            if missing(`_hminus') | `_hminus' <= 0 {
                di as error "h0:_cons missing or non-positive in e(b)"
                exit 111
            }
            local _tau = (1-`t0')/(1-`t1')
            local _x   = `_tau'^`elasticity'
            local _rho = ln(`_x')
            local _sl  = (`zlo' + `zl')/2 - `zstar'
            local _sr  = (`zh'  + `zhi')/2 - `zstar'
            if `islog_val' {
                local _sr0 = `_sr' + `_rho'
                local _hr0 = `_hplus'
            }
            else {
                local _sr0 = (`_x' - 1)*`zstar' + `_x'*`_sr'
                local _hr0 = `_hplus'/`_x'
            }
            // slope / intercept of the implicit two-point counterfactual
            // line, in running-variable coordinates
            local _den = `_sr0' - `_sl'
            local _m0 = 0
            if abs(`_den') >= 1e-12 {
                local _m = (`_hr0' - `_hminus')/`_den'
                local _a = `_hminus' - `_m'*`_sl'
                local _m0 = cond(`_a' > 0, `_m'/`_a', 0)
            }
            matrix `bcoef' = (1, `_m0')
        }
        else {
            // Polynomial estimator: pull all h0 coefficients from e(b) and
            // recentre them on the cutoff (in the running variable -- log
            // earnings when e(log)==1, since polbunch's `log' means the
            // variable is already logged).
            local _K      = e(polynomial)
            local _cest   = e(cutoff_est)
            local _xscale = e(xscale)
            matrix `craw' = J(1, `_K'+1, 0)
            matrix `craw'[1,1] = _b[h0:_cons]
            local _zterm "c.`_zname'"
            forvalues _k = 1/`_K' {
                if `_k' > 1 local _zterm "`_zterm'#c.`_zname'"
                capture local _bk = _b[h0:`_zterm']
                if !_rc & !missing(`_bk') matrix `craw'[1, `_k'+1] = `_bk'
            }
            mata: st_matrix("`bcoef'", pbx_ecoef_transform( ///
                st_matrix("`craw'"), `_cest', `_xscale'))
        }
    }
    else if `ncore' == 9 & ("`relslope'" != "" | "`h0poly'" != "") {
        // ================= EXPLICIT MODE =================
        local islog_val = ("`log'" != "")
        if "`bw'" == "" local bw = 1
        if "`h0poly'" != "" {
            local _nb : word count `h0poly'
            matrix `bcoef' = J(1, `_nb', 0)
            local _jb = 0
            foreach _v of local h0poly {
                local ++_jb
                matrix `bcoef'[1,`_jb'] = `_v'
            }
        }
        else {
            matrix `bcoef' = (1, `relslope')
        }
    }
    else {
        local missing_opts
        foreach v in estimator zstar t0 t1 elasticity zlo zhi zl zh {
            if "``v''" == "" local missing_opts `missing_opts' `v'()
        }
        di as error "Specify all 9 primary options plus relslope()/h0poly(), or none to use e()"
        di as error "Missing: `missing_opts' relslope()|h0poly()"
        exit 198
    }

    if !inlist(`estimator', 0, 1, 2, 3, 4) {
        di as error "polbunchbias: estimator must be 0, 1, 2, 3 or 4"
        exit 198
    }

    local _K1 = colsof(`bcoef')
    local _K  = `_K1' - 1

    // ---- resolve the two toggle axes ----
    if "`constant'" != ""         local useconstant = 1
    else if "`exact'" != ""       local useconstant = 0
    else if "`_e_constant'" != "" local useconstant = `_e_constant'
    else                          local useconstant = inlist(`estimator',1,2)

    if "`poolmass'" != ""         local nosplit = 1
    else if "`splitmass'" != ""   local nosplit = 0
    else if "`_e_nosplit'" != ""  local nosplit = `_e_nosplit'
    else                          local nosplit = inlist(`estimator',1,2)
    if `estimator' == 1 local nosplit = 1

    if inlist(`estimator',0,3) & !`useconstant' & !`nosplit' {
        di as error "polbunchbias: estimator(`estimator') is consistent under exact integration and a split mass calc -- zero bias by construction. Specify constant, poolmass, or both; or use estimator(1), estimator(2), or estimator(4)."
        exit 198
    }

    // Saez (estimator 4) with a constant-density inversion uses a FLAT
    // counterfactual both to form the excess mass and to invert it into a
    // response -- the two stages agree, so the internal-consistency bias
    // this command reports is zero by construction (same situation as
    // estimators 0/3 above).  This is NOT a claim that Saez/constant is
    // unbiased: it is biased whenever the true counterfactual slopes over
    // the excluded window.  estimator(4) + exact keeps a non-trivial
    // internal check (the excess mass is built from a step counterfactual
    // but inverted against the implied line); a comparable bias across
    // estimators needs a common counterfactual, not this per-estimator one.
    if `estimator' == 4 & `useconstant' {
        di as error "polbunchbias: estimator(4) with constant has no internal-consistency bias -- Saez"
        di as error "             with a constant-density inversion applies a flat counterfactual in both"
        di as error "             the mass and the inversion step, so the bias is zero by construction."
        di as error "             This is NOT a statement that the Saez/constant estimator is unbiased"
        di as error "             (it is biased whenever the counterfactual slopes across the window)."
        di as error "             Use estimator(4) with exact for the step-vs-line internal check, or"
        di as error "             compare estimators against a common counterfactual density."
        exit 198
    }

    if "`bw'" == "" local bw = 1
    local _doiter = ("`iterate'" != "")
    local _allowneg = ("`allownegative'" != "")
    // h0>=0 is a check on the entire fitted FUNCTION shape, not a single
    // theoretically-bounded scalar like the elasticity -- a polynomial
    // approximation to a noisily-ESTIMATED density can legitimately dip
    // slightly negative in a low-density region without the true h0 being
    // invalid (the usual artifact of any smooth/polynomial density
    // approximation), and an MC validation against a known truth found
    // this check can erase the entire benefit of promotion, not just
    // filter out bad candidates.  Opt-in (h0check), not opt-out.
    local _checkh0  = ("`h0check'" != "")

    // w-space working variables: bcoef is expressed in units where 1 raw
    // z-unit (or, in log mode, 1 unit of log-earnings) equals 1/_xscale.
    // zstarw is the cutoff's ECONOMIC magnitude in those units (used for
    // the Delta-to-level-shift anchor and the relative-slope prefactor);
    // its ABSOLUTE POSITION is implicitly 0 since _lo_w.._hi_w are already
    // offsets from it.  _xscale=1 (h0poly/relslope/Saez/explicit mode)
    // makes all of this identical to the old raw-coordinate locals.
    local _zstarw = `zstar'/`_xscale'
    local _lo_w   = (`zlo'-`zstar')/`_xscale'
    local _hi_w   = (`zhi'-`zstar')/`_xscale'
    local _L_w    = (`zl' -`zstar')/`_xscale'
    local _H_w    = (`zh' -`zstar')/`_xscale'
    local _bw_w   = `bw'/`_xscale'

    // ---- run (optionally iterated), with polynomial-order fallback ------
    // The fitting-stage bias (estimators 1, 2) inverts a degree-K monomial
    // Gram on the fitting window.  On a wide / off-centre window at a high
    // order that design can be near-collinear enough that pbx_bias_core
    // returns a missing bias (its own ill-conditioning / plausibility
    // guards).  When that happens, retry with the best degree-(K-1) L2
    // approximation of the counterfactual over the same window, stepping
    // down until the bias is well defined.  The h0 shape is essentially
    // unchanged whenever the dropped high-order terms were noise; the note
    // below reports the order actually used and the relative L2 change.
    local _Kfull = `_K'
    local _Kused = `_Kfull'
    local _bias_l2err = 0
    tempname _bcuse
    matrix `_bcuse' = `bcoef'
    forvalues _kk = `_Kfull'(-1)1 {
        if `_kk' < `_Kfull' & "`nofallback'" == "" {
            mata: st_matrix("`_bcuse'", pbx_polreduce(st_matrix("`bcoef'"), ///
                `_lo_w', `_L_w', `_H_w', `_hi_w', `_kk'))
            mata: st_numscalar("_pbx_l2", pbx_pol_l2rel(st_matrix("`bcoef'"), ///
                st_matrix("`_bcuse'"), `_lo_w', `_L_w', `_H_w', `_hi_w'))
            local _bias_l2err = _pbx_l2
        }
        capture matrix drop _pbxs_core _pbxs_bbeta _pbxs_bcur
        capture scalar drop _pbxs_ecur _pbxs_iter _pbxs_conv _pbxs_diverged
        mata: pbx_bias_solve_ws(`estimator', `bmodel', `islog_val', `useconstant', ///
            `nosplit', `_zstarw', `t0', `t1', st_matrix("`_bcuse'"), `elasticity', ///
            `_lo_w', `_hi_w', `_L_w', `_H_w', `_bw_w', `_xscale', `_doiter', ///
            `tolerance', `maxiter', `underrelax', `_allowneg', `_checkh0', "_pbxs")
        local _Kused = `_kk'
        if !inlist(`estimator',1,2)          continue, break
        if "`nofallback'" != ""              continue, break
        // accept this order only if it yields a usable bias: the fitting
        // design inverted (bias_h, col 20) AND the response inversion
        // produced a finite elasticity bias (col 24).  A degree that clears
        // the first but not the second is still degenerate -- keep stepping
        // down; if none qualifies the "not identified" path below fires.
        if !missing(el(_pbxs_core,1,20)) & !missing(el(_pbxs_core,1,24)) ///
                                            continue, break
    }
    local _bdiv0 = _pbxs_diverged
    capture scalar drop _pbx_l2

    tempname out bbeta bcur
    matrix `out'   = _pbxs_core
    matrix `bbeta' = _pbxs_bbeta
    matrix `bcur'  = _pbxs_bcur

    // ---- diverged iterate: is there a LOWER order where it converges? ----
    // The un-iterated bias above plugs the fitted (biased) h0/elasticity in
    // as if they were truth -- a leading-order approximation, not a fix for
    // the circularity iterate exists to solve.  A degree-K' < K fit that
    // actually reaches the self-consistent fixed point directly addresses
    // that circularity, at the cost of truncating h0 to a lower degree (the
    // BEST L2 approximation of the fitted h0, so only as much as needed).
    // A hand-fitted or digitized h0 is itself noisy at high order (its
    // highest-order terms are typically the least precisely estimated, or
    // outright insignificant) and polynomial extrapolation -- exactly what
    // the bias integral does across the excluded window -- amplifies that
    // noise; non-convergence at the full order is a symptom of the same
    // thing.  So by DEFAULT the converged lower-order fit is PROMOTED to the
    // primary reported bias, trading a small deliberate truncation for a lot
    // less extrapolation variance; the un-iterated full-order plug-in number
    // is kept, not discarded, as bias_uniter_*.  `nopromote' keeps the OLD
    // default instead (full-order plug-in primary, converged fit as
    // bias_iter_*); `nofallback' disables this search entirely.
    local _iterK = .
    local _iterl2 = .
    if `_doiter' & inlist(`estimator',1,2) & "`nofallback'" == "" & ///
            "`_bdiv0'" == "1" & `_Kused' > 1 {
        tempname _bciter _iterout _iterbbeta _iterbcur
        forvalues _kk2 = `=`_Kused'-1'(-1)1 {
            mata: st_matrix("`_bciter'", pbx_polreduce(st_matrix("`bcoef'"), ///
                `_lo_w', `_L_w', `_H_w', `_hi_w', `_kk2'))
            capture matrix drop _pbxi_core _pbxi_bbeta _pbxi_bcur
            capture scalar drop _pbxi_ecur _pbxi_iter _pbxi_conv _pbxi_diverged
            mata: pbx_bias_solve_ws(`estimator', `bmodel', `islog_val', `useconstant', ///
                `nosplit', `_zstarw', `t0', `t1', st_matrix("`_bciter'"), `elasticity', ///
                `_lo_w', `_hi_w', `_L_w', `_H_w', `_bw_w', `_xscale', 1, ///
                `tolerance', `maxiter', `underrelax', `_allowneg', `_checkh0', "_pbxi")
            if _pbxi_conv == 1 & _pbxi_diverged == 0 & ///
                    !missing(el(_pbxi_core,1,20)) & !missing(el(_pbxi_core,1,24)) {
                local _iterK = `_kk2'
                mata: st_numscalar("_pbx_l2b", pbx_pol_l2rel(st_matrix("`bcoef'"), ///
                    st_matrix("`_bciter'"), `_lo_w', `_L_w', `_H_w', `_hi_w'))
                local _iterl2 = _pbx_l2b
                matrix `_iterout'   = _pbxi_core
                matrix `_iterbbeta' = _pbxi_bbeta
                matrix `_iterbcur'  = _pbxi_bcur
                local _iterconv = _pbxi_conv
                local _iteriter = _pbxi_iter
                local _iterecur = _pbxi_ecur
                continue, break
            }
        }
        capture scalar drop _pbx_l2b
        capture matrix drop _pbxi_core _pbxi_bbeta _pbxi_bcur
        capture scalar drop _pbxi_ecur _pbxi_iter _pbxi_conv _pbxi_diverged
    }

    // ---- promote the converged lower-order fit to the primary result -----
    // GUARDRAIL: a converged fixed point is only worth promoting if it is
    // still recognisably the SAME counterfactual.  When the search has to
    // step all the way down to collapse most of h0's curvature to reach
    // convergence (large L2 change), the "self-consistent" fit is really a
    // different, over-simplified model, not a refinement of the original --
    // promoting it can make the correction WORSE, not better (confirmed by
    // MC: estimator 2's convergence criterion is harder to satisfy than
    // estimator 1's and often only converges after collapsing to K=1,
    // L2~0.6-0.7, which is a materially worse bias correction than the
    // un-iterated full-order fit; estimator 1 typically converges after a
    // modest step, L2~0.15-0.20, which IS an improvement).  promotetol()
    // sets the cutoff; default 0.3.
    local _promoted = 0
    local _uKused = .
    if `_iterK' < . & "`nopromote'" == "" & `_iterl2' <= `promotetol' {
        tempname _uout _ubbeta _ubcur
        matrix `_uout'  = `out'
        matrix `_ubbeta' = `bbeta'
        matrix `_ubcur'  = `bcur'
        local _uKused = `_Kused'
        matrix `out'   = `_iterout'
        matrix `bbeta' = `_iterbbeta'
        matrix `bcur'  = `_iterbcur'
        local _Kused   = `_iterK'
        local _promoted = 1
    }

    // fitting estimators always carry a fitting-stage bias; a missing bias_h
    // (or elasticity bias) means the design was singular / ill-conditioned on
    // the fitted window at every order down to 1 (estimator 2: or the right
    // window [`zh',`zhi'] alone is too short to identify Chetty's delta)
    if inlist(`estimator',1,2) & (missing(el(`out',1,20)) | missing(el(`out',1,24))) {
        di as error "polbunchbias: the counterfactual polynomial is not identified on the fitted window [`zlo',`zl'] u [`zh',`zhi'] at any order down to 1 -- window too short, or too collinear, for a bias estimate. Change the estimation window."
        // the fallback failed: report the full order, all-missing, so nothing
        // downstream shows a partial (degree-1) number
        local _Kused = `_Kfull'
        forvalues _wj = 20/26 {
            matrix `out'[1,`_wj'] = .
        }
        matrix `bbeta' = J(1, `_Kfull' + 1, .)
        matrix `bcur'  = J(1, `_Kfull' + 1, .)
    }
    else if inlist(`estimator',1,2) & cond(`_promoted',`_uKused',`_Kused') < `_Kfull' {
        // this is the CONDITIONING-driven reduction (before any promotion);
        // report it against the order that survived it, `_uKused' if promoted
        di as text "Note: the degree-`_Kfull' fitting design is too ill-conditioned on this window for a"
        di as text "      reliable bias; computed at degree " cond(`_promoted',`_uKused',`_Kused') ///
            " (best L2 fit to the degree-`_Kfull'"
        di as text "      counterfactual, relative change " %6.4f `_bias_l2err' ///
            cond(`_bias_l2err' < 0.02, " -- the same quantity).", " -- treat as indicative).")
    }

    // order actually used for the bias -> keep every downstream label,
    // betaname and the returned polynomial order consistent with it
    local _K  = `_Kused'
    local _K1 = `_Kused' + 1

    // pbx_bias_core_ws echoes columns 13-16 as offsets from the cutoff (its
    // own native w-space convention); r(zlo)/r(zhi)/r(zL)/r(zH) have always
    // meant the ABSOLUTE window bounds (matching explicit mode's own
    // arguments) -- restore that convention here, same patch the raw-
    // coordinate pbx_bias_core() wrapper applies for direct Mata callers.
    matrix `out'[1,13] = `zlo'
    matrix `out'[1,14] = `zhi'
    matrix `out'[1,15] = `zl'
    matrix `out'[1,16] = `zh'

    // "bmodel" is a vestigial always-0 column (the option was removed); kept
    // so the 26-column return layout and its consumers stay unchanged.
    local names estimator bmodel islog zstar t0 t1 tau lambda elasticity x rho Delta ///
        zlo zhi zL zH dL dR B bias_h bias_B bias_response bias_shift ///
        bias_elasticity bias_slope bias_lambda
    matrix colnames `out' = `names'

    // per-coefficient bias labels b0, b1, ...
    local betanames
    forvalues j = 0/`_K' {
        local betanames `betanames' b`j'
    }
    matrix colnames `bbeta' = `betanames'
    matrix colnames `bcur'  = `betanames'

    // display / r(b) matrix: b0..bK then the downstream estimands
    tempname rb tail
    matrix `tail' = ( el(`out',1,26), el(`out',1,21), el(`out',1,22), ///
        el(`out',1,23), el(`out',1,24) )
    matrix `rb' = `bbeta' , `tail'
    matrix colnames `rb' = `betanames' relative_slope number_bunchers ///
        marginal_response shift elasticity

    if `_promoted'          local _ordertxt "polynomial order `_Kused' (converged iterate; degree-`_uKused' plug-in did not converge)"
    else if `_Kused' < `_Kfull' local _ordertxt "polynomial order `_Kused' (fallback from `_Kfull')"
    else                    local _ordertxt "polynomial order `_K'"
    di as text _newline "Polynomial bunching bias estimates: estimator `estimator', `_ordertxt'"
    di as text "{hline 55}"
    di as text %25s "Estimand" "  " %12s "Bias"
    di as text "{hline 55}"
    local _dn : colnames `rb'
    local _dk = colsof(`rb')
    forvalues j = 1/`_dk' {
        local _dnm : word `j' of `_dn'
        di as text %25s "`_dnm'" "  " as result %12.6g `rb'[1,`j']
    }
    di as text "{hline 55}"
    if `estimator' == 4 local _reftxt "Saez two-point counterfactual line"
    else                local _reftxt "fitted counterfactual (degree `_K')"
    di as text "Note: bias is measured against this estimator's own"
    di as text "      `_reftxt'."
    di as text "      Each estimator assumes a different counterfactual, so"
    di as text "      these figures are NOT directly comparable across"
    di as text "      estimators.  For a cross-estimator comparison, evaluate"
    di as text "      every estimator against one common counterfactual density."

    // ---- returns ----
    forvalues j = 1/26 {
        local nm : word `j' of `names'
        tempname s`j'
        scalar `s`j'' = el(`out',1,`j')
        return scalar `nm' = `s`j''
    }
    return matrix b = `rb'
    return matrix bias_beta = `bbeta'
    return matrix corrected_beta = `bcur'
    return scalar polynomial = `_K'
    return scalar bias_polynomial = `_Kused'
    return scalar bias_polyfull   = `_Kfull'
    return scalar bias_l2err      = `_bias_l2err'
    return scalar constant = `useconstant'
    return scalar nosplit = `nosplit'
    return scalar input_elasticity = `elasticity'
    if `_promoted' {
        return scalar corrected_elasticity = `_iterecur'
        return scalar iterations = `_iteriter'
        return scalar converged = `_iterconv'
        return scalar iterate_diverged = 0
    }
    else {
        return scalar corrected_elasticity = _pbxs_ecur
        return scalar iterations = _pbxs_iter
        return scalar converged = _pbxs_conv
        return scalar iterate_diverged = _pbxs_diverged
    }

    if `_doiter' & _pbxs_diverged == 1 & !`_promoted' {
        di as error "Warning: the iterate self-consistency loop did not converge for this cell;"
        di as error "         reporting the un-iterated (non-self-consistent) bias. r(iterate_diverged)=1."
    }

    if `_promoted' {
        local _u_be = el(`_uout',1,24)
        di as text "Note: iterate did not converge at degree `_uKused' (the plug-in bias there was"
        di as text "      " %9.4f `_u_be' "); PROMOTED to the converged, self-consistent fit at degree"
        di as text "      `_Kused' shown above (L2 change " %6.4f `_iterl2' " from the degree-`_uKused'"
        di as text "      counterfactual).  The un-iterated degree-`_uKused' plug-in bias is kept, not"
        di as text "      discarded -- see r(bias_uniter_*).  Specify {cmd:nopromote} to report the"
        di as text "      plug-in number as primary instead (the pre-2026-09-04 default)."
        // bias_l2err is redefined here to the TOTAL change from the originally
        // requested degree-`_Kfull' counterfactual to the final (promoted,
        // converged) order -- pbx_pol_l2rel was evaluated against `bcoef',
        // not the intermediate conditioning-fallback reduction
        return scalar bias_l2err = `_iterl2'
        return scalar bias_uniter_polynomial  = `_uKused'
        return scalar bias_uniter_h           = el(`_uout',1,20)
        return scalar bias_uniter_B           = el(`_uout',1,21)
        return scalar bias_uniter_response    = el(`_uout',1,22)
        return scalar bias_uniter_shift       = el(`_uout',1,23)
        return scalar bias_uniter_elasticity  = `_u_be'
    }
    else if `_iterK' < . {
        // nopromote: report the converged fit as a secondary number instead
        local _iter_bh = el(`_iterout',1,20)
        local _iter_bB = el(`_iterout',1,21)
        local _iter_br = el(`_iterout',1,22)
        local _iter_bs = el(`_iterout',1,23)
        local _iter_be = el(`_iterout',1,24)
        di as text "Note: iterate did not converge at degree `_Kused', but DOES converge at degree"
        di as text "      `_iterK' (a genuinely self-consistent fixed point; L2 change " %6.4f `_iterl2' "):"
        di as text "      elasticity bias " %9.4f `_iter_be' " there, vs " %9.4f el(`out',1,24) ///
            " un-iterated at degree `_Kused' above."
        if "`nopromote'" != "" {
            di as text "      Not promoted: nopromote was specified."
        }
        else {
            di as text "      Not promoted: the L2 change (" %6.4f `_iterl2' ") exceeds promotetol(" ///
                %4.2f `promotetol' ") -- degree `_iterK' collapses too much of h0's curvature to"
            di as text "      trust as the primary estimate; raise promotetol() to accept it anyway."
        }
        return scalar bias_iter_polynomial = `_iterK'
        return scalar bias_iter_l2err      = `_iterl2'
        return scalar bias_iter_h          = `_iter_bh'
        return scalar bias_iter_B          = `_iter_bB'
        return scalar bias_iter_response   = `_iter_br'
        return scalar bias_iter_shift      = `_iter_bs'
        return scalar bias_iter_elasticity = `_iter_be'
    }

    capture matrix drop _pbxs_core _pbxs_bbeta _pbxs_bcur
    capture scalar drop _pbxs_ecur _pbxs_iter _pbxs_conv _pbxs_diverged
end

* ====================================================================
* Polynomial-algebra layer
*
* A polynomial is a real ROW vector b = (b0, b1, ..., bK), where b[k+1]
* is the coefficient on s^k.  Degree K = cols(b) - 1.  The running
* variable s is always centered at the cutoff: s = z - zstar.
*
* Naming: pbx_*() are the public helpers; pbx__*() are test-only.
* ====================================================================

capture mata: mata drop pbx_deg()
capture mata: mata drop pbx_polval()
capture mata: mata drop pbx_polvals()
capture mata: mata drop pbx_polderiv()
capture mata: mata drop pbx_polminval()
capture mata: mata drop pbx_polantideriv()
capture mata: mata drop pbx_polint()
capture mata: mata drop pbx_polmoment()
capture mata: mata drop pbx_polgram()
capture mata: mata drop pbx_illcond()
capture mata: mata drop pbx_monoscale()
capture mata: mata drop pbx_choose()
capture mata: mata drop pbx_polshift()
capture mata: mata drop pbx_polscale()
capture mata: mata drop pbx_polaffine()
capture mata: mata drop pbx_poladd()
capture mata: mata drop pbx_polsub()
capture mata: mata drop pbx_polmul()
capture mata: mata drop pbx_polcompose()
capture mata: mata drop pbx_poltrim()
capture mata: mata drop pbx_ecoef_transform()
capture mata: mata drop pbx_lambda()
capture mata: mata drop pbx_solvequad()
capture mata: mata drop pbx_polroot()
capture mata: mata drop pbx_lambda_report()
capture mata: mata drop pbx_bias_core()
capture mata: mata drop pbx__chk()
capture mata: mata drop pbx__chkv()
capture mata: mata drop pbx_selftest()

mata:

// -------------------------------------------------------------------
// degree of b
// -------------------------------------------------------------------
real scalar pbx_deg(real rowvector b)
{
    return(cols(b) - 1)
}

// -------------------------------------------------------------------
// Horner evaluation at a scalar
// -------------------------------------------------------------------
real scalar pbx_polvals(real rowvector b, real scalar s)
{
    real scalar k, K, acc
    K   = cols(b) - 1
    acc = b[K+1]
    for (k = K-1; k >= 0; k--) acc = acc*s + b[k+1]
    return(acc)
}

// -------------------------------------------------------------------
// Horner evaluation, elementwise over a real matrix of s values
// (returns a matrix of the same shape)
// -------------------------------------------------------------------
real matrix pbx_polval(real rowvector b, real matrix s)
{
    real scalar k, K
    real matrix acc
    K   = cols(b) - 1
    acc = J(rows(s), cols(s), b[K+1])
    for (k = K-1; k >= 0; k--) acc = acc :* s :+ b[k+1]
    return(acc)
}

// -------------------------------------------------------------------
// coefficients of the derivative p'(s)   (degree K-1; (0) if K==0)
// -------------------------------------------------------------------
real rowvector pbx_polderiv(real rowvector b)
{
    real scalar K, k
    real rowvector d
    K = cols(b) - 1
    if (K == 0) return((0))
    d = J(1, K, 0)
    for (k = 1; k <= K; k++) d[k] = k*b[k+1]
    return(d)
}

// -------------------------------------------------------------------
// minimum value of polynomial b over [lo,hi]: the two endpoints and
// every REAL root of b' strictly inside the interval (a degree>=2 poly
// can dip below its endpoint values in between).  Mata's built-in
// polyroots() shares pbx_polvals' ascending (b[k+1] = coeff of s^k)
// convention, so b and its derivative feed it directly -- no separate
// root-finder needed.  Used to check h0 >= 0 over the analysis window.
// -------------------------------------------------------------------
real scalar pbx_polminval(real rowvector b, real scalar lo, real scalar hi)
{
    real rowvector d
    real matrix rts
    real scalar mv, i, r, v
    mv = min((pbx_polvals(b, lo), pbx_polvals(b, hi)))
    if (cols(b) - 1 >= 2) {
        d = pbx_polderiv(b)
        if (cols(d) >= 2) {
            rts = polyroots(d)
            for (i = 1; i <= cols(rts); i++) {
                if (Im(rts[i]) == 0) {
                    r = Re(rts[i])
                    if (r > lo & r < hi) {
                        v = pbx_polvals(b, r)
                        if (v < mv) mv = v
                    }
                }
            }
        }
    }
    return(mv)
}

// -------------------------------------------------------------------
// coefficients of the antiderivative with zero integration constant:
//   B(s) = sum_k b_k s^(k+1)/(k+1)   ->   (0, b0, b1/2, b2/3, ...)
// (degree K+1)
// -------------------------------------------------------------------
real rowvector pbx_polantideriv(real rowvector b)
{
    real scalar K, k
    real rowvector B
    K = cols(b) - 1
    B = J(1, K+2, 0)
    for (k = 0; k <= K; k++) B[k+2] = b[k+1] / (k+1)
    return(B)
}

// -------------------------------------------------------------------
// definite integral  int_lo^hi poly_b(s) ds
// -------------------------------------------------------------------
real scalar pbx_polint(real rowvector b, real scalar lo, real scalar hi)
{
    real scalar K, k, acc
    K   = cols(b) - 1
    acc = 0
    for (k = 0; k <= K; k++)
        acc = acc + b[k+1]*(hi^(k+1) - lo^(k+1))/(k+1)
    return(acc)
}

// -------------------------------------------------------------------
// moment vector: int_lo^hi P_Kout(s) * poly_b(s) ds
//   returns a (Kout+1) COLUMN vector; entry i (i = 0..Kout) is
//     sum_j b_j (hi^{i+j+1} - lo^{i+j+1}) / (i+j+1)
// (Kout is the OUTPUT polynomial degree -- may differ from deg(b))
// -------------------------------------------------------------------
real colvector pbx_polmoment(real rowvector b, real scalar lo,
                             real scalar hi, real scalar Kout)
{
    real scalar i, j, Kb, acc
    real colvector v
    Kb = cols(b) - 1
    v  = J(Kout+1, 1, 0)
    for (i = 0; i <= Kout; i++) {
        acc = 0
        for (j = 0; j <= Kb; j++)
            acc = acc + b[j+1]*(hi^(i+j+1) - lo^(i+j+1))/(i+j+1)
        v[i+1] = acc
    }
    return(v)
}

// -------------------------------------------------------------------
// Gram matrix: int_lo^hi P_K(s) P_K(s)' ds
//   (K+1)x(K+1), [i,j] = (hi^{i+j+1} - lo^{i+j+1}) / (i+j+1)
// -------------------------------------------------------------------
real matrix pbx_polgram(real scalar lo, real scalar hi, real scalar K)
{
    real scalar i, j
    real matrix G
    G = J(K+1, K+1, 0)
    for (i = 0; i <= K; i++)
        for (j = 0; j <= K; j++)
            G[i+1,j+1] = (hi^(i+j+1) - lo^(i+j+1))/(i+j+1)
    return(G)
}

// -------------------------------------------------------------------
// is a symmetric PSD design matrix too ill-conditioned to invert
// safely at this polynomial order?  (fitting window too short for K)
// -------------------------------------------------------------------
real scalar pbx_illcond(real matrix M)
{
    real rowvector ev
    real scalar lo, hi
    ev = symeigenvalues(M)
    lo = min(ev)
    hi = max(ev)
    if (hi <= 0) return(1)
    return(lo <= hi*1e-12)
}

// -------------------------------------------------------------------
// diagonal preconditioner for a degree-K monomial (Gram / moment)
// design on the centred running variable.  Returns the (K+1) column
// vector d with d[k+1] = sc^{-k}, sc = the largest |window endpoint|
// (floored at 1).  Then (d d') :* (int s^{i+j} ds) has O(1) entries,
// so a wide fitting window (|s| ~ 50, e.g. a cutoff far from 0 or a
// long excluded region) no longer blows the condition number and
// trips pbx_illcond().  Rescaling the linear system by d and its
// solution back by d is an exact change of basis.
// -------------------------------------------------------------------
real colvector pbx_monoscale(real scalar lo, real scalar hi,
    real scalar L, real scalar H, real scalar K)
{
    real scalar sc
    sc = max((abs(lo), abs(hi), abs(L), abs(H), 1))
    return(sc :^ (-(0::K)))
}

// -------------------------------------------------------------------
// best degree-Kt L2 approximation of a degree-K polynomial b over the
// fitting window [lo,L] u [H,hi] (offsets from the cutoff).  Used by
// the polynomial-order fallback: when the degree-K fitting design is
// too ill-conditioned for a reliable bias, retry with this reduced
// counterfactual.  Same diagonal preconditioning as the bias engine.
// -------------------------------------------------------------------
real rowvector pbx_polreduce(real rowvector b, real scalar lo,
    real scalar L, real scalar H, real scalar hi, real scalar Kt)
{
    real matrix G
    real colvector m, dsc, c
    G   = pbx_polgram(lo, L, Kt) + pbx_polgram(H, hi, Kt)
    m   = pbx_polmoment(b, lo, L, Kt) + pbx_polmoment(b, H, hi, Kt)
    dsc = pbx_monoscale(lo, hi, L, H, Kt)
    c   = dsc :* (invsym((dsc * dsc') :* G) * (dsc :* m))
    return(c')
}

// -------------------------------------------------------------------
// relative L2 distance ||bf - br|| / ||bf|| over [lo,L] u [H,hi].
// Reports how much the order reduction changed the counterfactual.
// -------------------------------------------------------------------
real scalar pbx_pol_l2rel(real rowvector bf, real rowvector br,
    real scalar lo, real scalar L, real scalar H, real scalar hi)
{
    real rowvector d, dd, ff
    real scalar num, den
    d   = pbx_polsub(bf, br)
    dd  = pbx_polmul(d, d)
    ff  = pbx_polmul(bf, bf)
    num = pbx_polint(dd, lo, L) + pbx_polint(dd, H, hi)
    den = pbx_polint(ff, lo, L) + pbx_polint(ff, H, hi)
    if (den <= 0) return(.)
    return(sqrt(num/den))
}

// -------------------------------------------------------------------
// exact binomial coefficient C(n,k) for small non-negative integers
// -------------------------------------------------------------------
real scalar pbx_choose(real scalar n, real scalar k)
{
    real scalar r, i
    if (k < 0 | k > n) return(0)
    r = 1
    for (i = 1; i <= k; i++) r = r*(n - k + i)/i
    return(r)
}

// -------------------------------------------------------------------
// coefficients of p(s + c)   (same degree; Taylor coefficients of p at c)
//
// Computed by REPEATED SYNTHETIC DIVISION ("Horner shift"), not the
// direct binomial-coefficient sum p_j = sum_{k>=j} b_k C(k,j) c^(k-j).
// Both give the same answer in exact arithmetic and cost O(K^2) either
// way, but the direct sum forms large intermediate terms (binomial
// coefficients times powers of c) that can nearly cancel -- exactly the
// classic ill-conditioning of a monomial-basis recentring at high degree,
// and the actual source of precision loss behind polbunch's e()-mode /
// iterate instability at high polynomial order (rescaling afterwards, as
// pbx_monoscale does elsewhere, cannot recover precision already lost
// here).  Synthetic division instead peels off one Taylor coefficient at
// a time as the remainder of dividing by (s - c), each step a simple
// a + c*b update with no large cancelling sum -- the standard numerically
// stable algorithm for this operation.
// -------------------------------------------------------------------
real rowvector pbx_polshift(real rowvector b, real scalar c)
{
    real scalar K, i, k, n
    real rowvector a, bq, out
    K = cols(b) - 1
    if (K == 0) return(b)
    out = J(1, K+1, 0)
    a   = b
    n   = K + 1
    for (i = 0; i <= K; i++) {
        if (n == 1) {
            out[i+1] = a[1]
            continue
        }
        bq = J(1, n-1, 0)
        bq[n-1] = a[n]
        for (k = n-2; k >= 1; k--) bq[k] = a[k+1] + c*bq[k+1]
        out[i+1] = a[1] + c*bq[1]
        a = bq
        n = n - 1
    }
    return(out)
}

// -------------------------------------------------------------------
// coefficients of p(A*s)   ->   new b_k = b_k * A^k
// -------------------------------------------------------------------
real rowvector pbx_polscale(real rowvector b, real scalar A)
{
    real scalar K, k
    real rowvector out
    K   = cols(b) - 1
    out = b
    for (k = 0; k <= K; k++) out[k+1] = out[k+1]*A^k
    return(out)
}

// -------------------------------------------------------------------
// coefficients of  A * p(A*s + c)
//   level-case h1 relocation:  h1(s) = x * h0(x*s + r)
//   ->  pbx_polaffine(h0, x, r)
// -------------------------------------------------------------------
real rowvector pbx_polaffine(real rowvector b, real scalar A, real scalar c)
{
    return(A :* pbx_polscale(pbx_polshift(b, c), A))
}

// -------------------------------------------------------------------
// add / subtract polynomials of possibly different degree
// -------------------------------------------------------------------
real rowvector pbx_poladd(real rowvector a, real rowvector b)
{
    real scalar na, nb, n, k
    real rowvector out
    na = cols(a); nb = cols(b)
    n  = max((na, nb))
    out = J(1, n, 0)
    for (k = 1; k <= na; k++) out[k] = out[k] + a[k]
    for (k = 1; k <= nb; k++) out[k] = out[k] + b[k]
    return(out)
}

real rowvector pbx_polsub(real rowvector a, real rowvector b)
{
    return(pbx_poladd(a, -b))
}

// -------------------------------------------------------------------
// polynomial product  (a*b)_n = sum_{i+j=n} a_i b_j   (degree da+db)
// -------------------------------------------------------------------
real rowvector pbx_polmul(real rowvector a, real rowvector b)
{
    real scalar na, nb, i, j
    real rowvector out
    na = cols(a) ; nb = cols(b)
    out = J(1, na + nb - 1, 0)
    for (i = 1; i <= na; i++)
        for (j = 1; j <= nb; j++)
            out[i+j-1] = out[i+j-1] + a[i]*b[j]
    return(out)
}

// -------------------------------------------------------------------
// composition  a(b(s))  truncated to degree Kout
// -------------------------------------------------------------------
real rowvector pbx_polcompose(real rowvector a, real rowvector b,
                              real scalar Kout)
{
    real scalar Ka, k, j
    real rowvector out, pw
    Ka  = cols(a) - 1
    out = J(1, Kout+1, 0)
    out[1] = a[1]                       // a_0 * b(s)^0
    pw = (1)                            // b(s)^0
    for (k = 1; k <= Ka; k++) {
        pw = pbx_polmul(pw, b)
        if (cols(pw) > Kout+1) pw = pw[|1 \ Kout+1|]
        for (j = 1; j <= cols(pw); j++) out[j] = out[j] + a[k+1]*pw[j]
    }
    return(out)
}

// -------------------------------------------------------------------
// trim trailing (near-)zero coefficients; keeps at least the constant
// -------------------------------------------------------------------
real rowvector pbx_poltrim(real rowvector b, real scalar tol)
{
    real scalar k
    k = cols(b)
    while (k > 1 & abs(b[k]) <= tol) k--
    return(b[1..k])
}

// -------------------------------------------------------------------
// NORMALISED-SPACE VERSION.  polbunch fits h0 as coefficients c of its
// own normalised coordinate u = (z - zmid)/xscale.  The production
// pbx_ecoef_transform() shifts to the cutoff (u -> u - cest) AND rescales
// out of u-space into raw z-units (s = z - z*) in the same step -- that
// rescale is exactly what stretches the coefficient vector across many
// orders of magnitude when xscale is far from 1 (b0 ~ 1, a high-degree
// b_K ~ xscale^-K), which is the root cause this file exists to remove.
//
// Here we do ONLY the shift (still numerically the delicate part, hence
// the synthetic-division pbx_polshift), staying in u-space: the returned
// polynomial is h0 as a function of w = u - cest = (z - z*)/xscale, i.e.
// distance from the cutoff MEASURED IN polbunch's OWN normalised units,
// never converted to raw z-units at all.  pbx_bias_core is written to
// consume w-space coefficients directly (see its header comment for the
// unit bookkeeping this implies for zstar, bw, and the outputs).
// -------------------------------------------------------------------
real rowvector pbx_ecoef_transform(real rowvector c, real scalar cest,
    real scalar xscale)
{
    real rowvector braw
    braw = pbx_polshift(c, cest)
    if (abs(braw[1]) > 1e-300) braw = braw :/ braw[1]
    return(braw)
}

// -------------------------------------------------------------------
// local relative slope  lambda = zstar * h0'(zstar)/h0(zstar)
// from centered coefficients (needs degree >= 1)
// -------------------------------------------------------------------
real scalar pbx_lambda(real rowvector b, real scalar zstar)
{
    if (cols(b) < 2)          return(.)
    if (abs(b[1]) < 1e-14)    return(.)
    return(zstar * b[2] / b[1])
}

// -------------------------------------------------------------------
// solve   0.5*mm*r^2 + aa*r = S   for r
// (mirrors pb_solve_quad() of polbunchbias.ado exactly, so the
//  degree-1 path reproduces the current command bit-for-bit)
// -------------------------------------------------------------------
real scalar pbx_solvequad(real scalar aa, real scalar mm, real scalar S)
{
    real scalar disc, r1, r2

    if (abs(mm) < 1e-12) {
        if (abs(aa) < 1e-12) return(.)
        return(S/aa)
    }
    disc = aa^2 + 2*mm*S
    if (disc < 0) return(.)

    r1 = (-aa + sqrt(disc))/mm
    r2 = (-aa - sqrt(disc))/mm

    if (r1 >= 0 & (r2 < 0 | r1 <= r2)) return(r1)
    if (r2 >= 0) return(r2)
    return(r1)
}

// -------------------------------------------------------------------
// solve   int_0^r poly_b(s) ds = rhs   for r
//
//   F(r) = sum_k b_k r^{k+1}/(k+1) - rhs
//
// deg(b) <= 1  -> analytic (pbx_solvequad).
// deg(b) >= 2  -> Newton from `guess` (or a constant-density guess),
//                 then a bisection sweep for the smallest positive
//                 root as a fallback.  Returns . if none is found.
// -------------------------------------------------------------------
real scalar pbx_polroot(real rowvector b, real scalar rhs, real scalar guess)
{
    real scalar K, r, fr, dfr, step, it
    real scalar cap, nseg, prev, cur, lo, hi, xa, xb, flo, fmid, mid, j
    real rowvector F, dF

    K = cols(b) - 1

    if (K == 0) {
        if (abs(b[1]) < 1e-14) return(.)
        return(rhs / b[1])
    }
    if (K == 1) {
        return(pbx_solvequad(b[1], b[2], rhs))
    }

    // F(r) and F'(r) as polynomials in r
    F    = pbx_polantideriv(b)      // length K+2, F[1] == 0
    F[1] = -rhs
    dF   = b

    // ---- initial guess: user guess if positive, else constant-density ----
    if (guess < . & guess > 0) {
        r = guess
    }
    else if (abs(b[1]) > 1e-14) {
        r = rhs / b[1]
    }
    else {
        r = 1
    }

    // ---- Newton ----
    for (it = 1; it <= 100; it++) {
        fr  = pbx_polvals(F, r)
        dfr = pbx_polvals(dF, r)
        if (abs(dfr) < 1e-14 | abs(r) > 1e12) {
            break
        }
        step = fr/dfr
        r = r - step
        if (abs(step) < 1e-13) {
            break
        }
    }
    if (!missing(r)) {
        if (r >= 0 & abs(pbx_polvals(F, r)) < 1e-7) {
            return(r)
        }
    }

    // ---- bisection sweep for the smallest positive root ----
    cap = 1
    if (guess < . & guess > 0) {
        cap = max((cap, 3*guess))
    }
    if (abs(b[1]) > 1e-14) {
        cap = max((cap, 3*abs(rhs/b[1])))
    }
    nseg = 400
    prev = pbx_polvals(F, 0)       // = -rhs
    lo   = 0
    for (it = 1; it <= nseg; it++) {
        hi  = cap*it/nseg
        cur = pbx_polvals(F, hi)
        if (prev == 0) {
            return(lo)
        }
        if (prev*cur < 0) {
            xa = lo
            xb = hi
            flo = prev
            for (j = 1; j <= 200; j++) {
                mid  = (xa + xb)/2
                fmid = pbx_polvals(F, mid)
                if (abs(fmid) < 1e-12 | (xb - xa) < 1e-14) {
                    return(mid)
                }
                if (flo*fmid < 0) {
                    xb = mid
                }
                else {
                    xa = mid
                    flo = fmid
                }
            }
            return((xa + xb)/2)
        }
        prev = cur
        lo   = hi
    }
    return(.)
}

// ===================================================================
// Test-only helpers
// ===================================================================

real scalar pbx__chk(string scalar desc, real scalar got,
                     real scalar want, real scalar tol)
{
    real scalar ok
    string scalar tag
    ok = (abs(got - want) <= tol)
    if (!ok & want != 0) ok = (abs(got/want - 1) <= tol)
    tag = "ok"
    if (!ok) tag = "<<<<< FAIL"
    printf("  %-44s got=%13.7g  want=%13.7g   %s\n", desc, got, want, tag)
    return(!ok)
}

real scalar pbx__chkv(string scalar desc, real matrix got, real matrix want,
                      real scalar tol)
{
    real scalar d, ok
    string scalar tag
    d  = max(abs(vec(got - want)))
    ok = (d <= tol)
    tag = "ok"
    if (!ok) tag = "<<<<< FAIL"
    printf("  %-44s maxabsdiff=%13.7g                %s\n", desc, d, tag)
    return(!ok)
}

// ===================================================================
// pbx_selftest(): exercises the polynomial-algebra layer against
// independently computed values.  Returns the number of failures.
// ===================================================================
real scalar pbx_selftest()
{
    real scalar f, x, r, rr, rhs
    real rowvector b, q, h0, h1
    real matrix S

    f = 0
    printf("\n{txt}pbx_selftest -- polynomial-algebra layer\n")
    printf("{hline 78}\n")

    // ---- evaluation: p(s) = 2 - 3s + s^2 ----
    b = (2, -3, 1)
    f = f + pbx__chk("polvals p(0)",  pbx_polvals(b, 0),  2, 1e-12)
    f = f + pbx__chk("polvals p(1)",  pbx_polvals(b, 1),  0, 1e-12)
    f = f + pbx__chk("polvals p(2)",  pbx_polvals(b, 2),  0, 1e-12)
    f = f + pbx__chk("polvals p(-1)", pbx_polvals(b, -1), 6, 1e-12)
    f = f + pbx__chkv("polval matrix vs polvals",
                      pbx_polval(b, (-1, 0, 1 \ 2, 3, 0.5)),
                      (6, 2, 0 \ 0, 2, 0.75), 1e-12)
    f = f + pbx__chk("deg", pbx_deg(b), 2, 0)

    // ---- derivative / antiderivative are inverses ----
    f = f + pbx__chkv("polderiv(2,-3,1)", pbx_polderiv(b), (-3, 2), 1e-12)
    f = f + pbx__chkv("polantideriv(2,-3,1)",
                      pbx_polantideriv(b), (0, 2, -1.5, 1/3), 1e-12)
    f = f + pbx__chkv("polderiv(polantideriv(b)) == b",
                      pbx_polderiv(pbx_polantideriv(b)), b, 1e-12)

    // ---- definite integrals ----
    f = f + pbx__chk("polint [0,1]  (2-3s+s^2)",
                     pbx_polint(b, 0, 1), 2 - 1.5 + 1/3, 1e-12)
    f = f + pbx__chk("polint [-1,2] (2-3s+s^2)",
                     pbx_polint(b, -1, 2), 4.5, 1e-12)
    f = f + pbx__chk("polint == polmoment(...,0)[1]",
                     pbx_polint(b, -0.7, 1.3),
                     pbx_polmoment(b, -0.7, 1.3, 0)[1], 1e-12)

    // ---- moment vector, unit polys reproduce the Gram matrix ----
    f = f + pbx__chkv("polmoment(1; [0,1]; K=2)",
                      pbx_polmoment((1), 0, 1, 2), (1 \ 0.5 \ 1/3), 1e-12)
    f = f + pbx__chkv("polmoment(s; [0,1]; K=1)",
                      pbx_polmoment((0, 1), 0, 1, 1), (0.5 \ 1/3), 1e-12)
    S = pbx_polgram(-0.3, 1.7, 3)
    f = f + pbx__chkv("gram col 0 == polmoment(e0)",
                      S[,1], pbx_polmoment((1), -0.3, 1.7, 3), 1e-10)
    f = f + pbx__chkv("gram col 2 == polmoment(e2)",
                      S[,3], pbx_polmoment((0, 0, 1), -0.3, 1.7, 3), 1e-10)
    f = f + pbx__chk("gram symmetric", max(abs(vec(S - S'))), 0, 1e-12)

    // ---- gram(.,.,1) equals the analytic 2x2 [ (h-l), (h^2-l^2)/2 ;
    //      . , (h^3-l^3)/3 ] used by the current command ----
    f = f + pbx__chkv("gram(l,h,1) analytic",
                      pbx_polgram(-0.3, 1.7, 1),
                      (1.7 - -0.3,             (1.7^2 - 0.09)/2 \
                       (1.7^2 - 0.09)/2,       (1.7^3 - -0.027)/3), 1e-12)

    // ---- argument shift:  p(s+c) ----
    //   p(s) = 2 - 3s + s^2 ;  p(s+1) = -s + s^2
    f = f + pbx__chkv("polshift(b,1)", pbx_polshift(b, 1), (0, -1, 1), 1e-12)
    q = pbx_polshift(b, 0.37)
    f = f + pbx__chk("polshift matches p(s+c) at s=0.9",
                     pbx_polvals(q, 0.9), pbx_polvals(b, 0.9 + 0.37), 1e-12)
    f = f + pbx__chk("polshift matches p(s+c) at s=-2.1",
                     pbx_polvals(q, -2.1), pbx_polvals(b, -2.1 + 0.37), 1e-12)

    // ---- argument scale:  p(A*s) ----
    f = f + pbx__chkv("polscale(b,2)", pbx_polscale(b, 2), (2, -6, 4), 1e-12)

    // ---- affine relocation  A*p(A*s + c)  (the level-case h1 map) ----
    //   linear check against the local-linear closed form:
    //   h0 = (a,m); h1 = ( x*a + x*m*r , x^2*m )
    h0 = (1, 0.5)
    x  = 1.2
    r  = 0.3
    f = f + pbx__chkv("polaffine linear vs closed form",
                      pbx_polaffine(h0, x, r), (x*1 + x*0.5*r, x^2*0.5), 1e-12)
    // general quadratic h0: verify A*p(A*s+c) pointwise
    h0 = (1, 0.3, -0.1)
    h1 = pbx_polaffine(h0, x, r)
    f = f + pbx__chk("polaffine pointwise s=0.8",
                     pbx_polvals(h1, 0.8), x*pbx_polvals(h0, x*0.8 + r), 1e-12)
    f = f + pbx__chk("polaffine pointwise s=-1.4",
                     pbx_polvals(h1, -1.4), x*pbx_polvals(h0, x*-1.4 + r), 1e-12)

    // ---- add / subtract ----
    f = f + pbx__chkv("poladd different degree",
                      pbx_poladd((1, 2), (0, 0, 3)), (1, 2, 3), 1e-12)
    f = f + pbx__chkv("polsub self == 0",
                      pbx_polsub((1, 2, 3), (1, 2, 3)), (0, 0, 0), 1e-12)

    // ---- polynomial product ----
    //   (2 - 3s + s^2)(1 + 2s) = 2 + s - 5s^2 + 2s^3
    f = f + pbx__chkv("polmul (2,-3,1)*(1,2)",
                      pbx_polmul((2, -3, 1), (1, 2)), (2, 1, -5, 2), 1e-12)
    f = f + pbx__chk("polmul matches pointwise product at s=1.3",
                     pbx_polvals(pbx_polmul((2,-3,1),(1,2)), 1.3),
                     pbx_polvals((2,-3,1),1.3)*pbx_polvals((1,2),1.3), 1e-12)
    f = f + pbx__chk("int (P'b)^2 == b' Gram b",
                     pbx_polint(pbx_polmul((1,0.4,-0.2),(1,0.4,-0.2)), -0.3, 1.1),
                     ((1,0.4,-0.2) * pbx_polgram(-0.3,1.1,2) * (1,0.4,-0.2)'),
                     1e-10)

    // ---- composition a(b(s)) truncated ----
    //   a = 1 + 2w - w^2 ,  b = s + 0.5 s^2  ->  a(b(s)) to deg 2:
    //   1 + 2(s+.5s^2) - (s+.5s^2)^2 = 1 + 2s + s^2 - s^2 - ... = 1 + 2s + 0*s^2
    f = f + pbx__chkv("polcompose deg2",
                      pbx_polcompose((1,2,-1), (0,1,0.5), 2), (1, 2, 0), 1e-12)
    f = f + pbx__chk("polcompose matches a(b(s)) at s=0.3",
                     pbx_polvals(pbx_polcompose((1,2,-1,0.5),(0,1,0.5,-0.1), 9), 0.3),
                     pbx_polvals((1,2,-1,0.5), pbx_polvals((0,1,0.5,-0.1),0.3)), 1e-12)

    // ---- relative slope ----
    f = f + pbx__chk("lambda z* m/a", pbx_lambda((2, 0.6, 9), 4),
                     4*0.6/2, 1e-12)

    // ---- root solve: int_0^r poly = rhs ----
    // K=0
    f = f + pbx__chk("polroot K=0", pbx_polroot((2), 6, .), 3, 1e-10)
    // K=1 (must match pbx_solvequad and hence pb_solve_quad)
    rhs = pbx_polint((1, 0.5), 0, 0.4)
    f = f + pbx__chk("polroot K=1 recovers r=0.4",
                     pbx_polroot((1, 0.5), rhs, .), 0.4, 1e-10)
    // K=2
    b   = (1, 0.3, -0.1)
    rhs = pbx_polint(b, 0, 0.5)
    f = f + pbx__chk("polroot K=2 recovers r=0.5 (guess)",
                     pbx_polroot(b, rhs, 0.5), 0.5, 1e-9)
    f = f + pbx__chk("polroot K=2 recovers r=0.5 (no guess)",
                     pbx_polroot(b, rhs, .), 0.5, 1e-9)
    // K=3
    b   = (1, 0.2, -0.05, 0.01)
    rhs = pbx_polint(b, 0, 0.7)
    rr  = pbx_polroot(b, rhs, .)
    f = f + pbx__chk("polroot K=3 recovers r=0.7", rr, 0.7, 1e-8)
    f = f + pbx__chk("polroot K=3 residual ~ 0",
                     pbx_polint(b, 0, rr) - rhs, 0, 1e-9)

    printf("{hline 78}\n")
    printf("{txt}pbx_selftest: %g check(s) FAILED\n\n", f)
    return(f)
}

end

* ====================================================================
* Bias core
*
* pbx_bias_core() computes the (exact, for a polynomial DGP) bias of the
* five estimators.  The counterfactual density h0 is a polynomial of ANY
* degree K in the centered running variable s, supplied as the row vector
* `bcoef` (bcoef[1] = height at the cutoff, bcoef[2] = running-coordinate
* slope, ...).  A degree-1 bcoef reproduces the earlier local-linear
* polbunchbias exactly.
*
* K = cols(bcoef) - 1 is BOTH the true h0 degree AND the polynomial order
* the estimator fits (they always match -- e()-mode reads e(polynomial)
* and the fitted h0 coefficients together).  Trailing zeros are therefore
* NOT inert: (1, m, 0, 0) models an estimator fitting a cubic to a truly
* linear h0, which carries more fitting-stage bias than (1, m).
*
* Return vector: 26 columns, then the (K+1) per-coefficient biases.
* The first 26 (SEE THE UNIT-BOOKKEEPING NOTE ABOVE pbx_bias_core --
* zstar/zlo/zhi/zL/zH/dL/dR here are RECONSTRUCTED z-unit echoes of the
* w-space working variables, not the original absolute inputs):
*   1  estimator     10 x            19 B             25 bias_slope
*   2  bmodel        11 rho          20 bias_h        26 bias_lambda
*   3  islog         12 Delta        21 bias_B
*   4  zstar         13 zlo          22 bias_response
*   5  t0            14 zhi          23 bias_shift
*   6  t1            15 zL           24 bias_elasticity
*   7  tau           16 zH
*   8  lambda        17 dL
*   9  elasticity    18 dR
* ====================================================================

capture mata: mata drop pbx_lambda_report()
capture mata: mata drop pbx_bias_core()
capture mata: mata drop pbx_bias_core_ws()
capture mata: mata drop pbx_bias_solve_ws()
capture mata: mata drop pbx_bias_scres_ws()
capture mata: mata drop pbx_bias_newton_ws()
capture mata: mata drop pbx_bias_e2_objQ()
capture mata: mata drop pbx_bias_e2_fixedpoint()
capture mata: mata drop pbx_bias_e2_solve()

mata:

// reported relative slope: (islog ? 1 : zstarw) * h0'(z*)/h0(z*).  Callers
// pass zstarw (= zstar/xscale, or zstar itself when xscale=1); the formula
// is unchanged -- see pbx_bias_core's header for why zstarw is the right
// value here (it is the "economic magnitude" role of the cutoff, distinct
// from its role as an absolute POSITION, which working in w-space already
// sets to 0 by construction).
// b's derivative b[2] is now taken w.r.t. w = s/xscale, so b[2]_w =
// b[2]_s * xscale for ANY b -- log mode's pref=1 (no zstar-style rescaling
// at all in the raw-coordinate code) must become pref=1/xscale to cancel
// that factor and recover the same coordinate-invariant number as before;
// level mode's zstarw already does the analogous cancellation (zstarw =
// zstar/xscale exactly compensates b[2]_w's extra factor of xscale).
real scalar pbx_lambda_report(real rowvector b, real scalar zstarw,
                              real scalar islog, real scalar xscale)
{
    real scalar pref
    if (cols(b) < 2)       return(0)
    if (abs(b[1]) < 1e-14) return(.)
    pref = zstarw
    if (islog) pref = 1/xscale
    return(pref * b[2] / b[1])
}

// ===================================================================
// NORMALISED-SPACE VERSION.
//
// bcoef is now a polynomial in w = (z - z*)/xscale -- polbunch's OWN
// normalised coordinate, re-centred at the cutoff but never rescaled back
// into raw z-units (see pbx_ecoef_transform).  Every argument that used to
// be a raw z-unit quantity is replaced by its w-space counterpart:
//
//   zstar (dual role in the old code)  ->  TWO separate things here:
//     - as an ABSOLUTE POSITION (window centring: lo=zlo-zstar etc.) this
//       is now simply 0 by construction, since w is already centred at
//       the cutoff -- the caller passes PRE-CENTRED bounds lo,hi,L,H
//       directly, so this role disappears from the function entirely.
//     - as an ECONOMIC MAGNITUDE (Delta -> level-shift anchor r=zstar*
//       Delta; the relative-slope prefactor) this is the argument named
//       zstarw below (= zstar/xscale).  Both roles happened to be the
//       SAME number in raw coordinates, which is exactly what made them
//       easy to conflate; they are not the same number here.
//   zlo,zhi,zL,zH  ->  lo,hi,L,H, already offsets from the cutoff, in
//       w-units (caller computes (zlo-zstar)/xscale etc.)
//   bw  ->  passed in ALREADY divided by xscale by the caller (bw_w).
//
// Dimensional bookkeeping (derived once, used throughout): bcoef's VALUES
// are coordinate-invariant (polbunch's h0 fit predicts a COUNT, not a
// density needing a Jacobian correction, so relabelling its argument axis
// does not rescale the function) -- bias_h, bias_B, bias_shift, bias_e and
// bias_lambda therefore come out correct with NO extra conversion, as long
// as zstarw and bw (already w-scaled) are used consistently in every
// formula that used to reference zstar/bw.  bias_response (a LEVEL shift)
// and bias_slope (a DERIVATIVE) do NOT come out unit-invariant -- they are
// computed in w-units throughout and converted back to z-units in a single
// place at the very end (bias_resp *xscale, bias_slope /xscale,
// biasbeta rescaled per-coefficient via pbx_polscale(.,1/xscale)), exactly
// mirroring what pbx_ecoef_transform used to do to the INPUT polynomial.
// The zstar/zlo/zhi/zL/zH/dL/dR echo columns are reconstructed from the
// w-space workings purely for human-readable display; nothing downstream
// consumes them.
// ===================================================================

// -------------------------------------------------------------------
// BACKWARD-COMPATIBILITY WRAPPER -- pre-2.2.0 interface: raw z-unit
// zstar and zlo/zhi/zL/zH (not yet centred on the cutoff), no xscale
// (equivalent to xscale=1, i.e. bcoef already in raw z-units).  This is
// exactly explicit/standalone mode's calling convention, so xscale=1 is
// the right value here, not a placeholder; direct Mata callers (the
// regression fixtures) use this signature and need no changes.  The .ado
// layer itself calls pbx_bias_core_ws/pbx_bias_solve_ws directly with the
// real (possibly !=1) xscale.
// -------------------------------------------------------------------
real rowvector pbx_bias_core(
    real scalar estimator, real scalar bmodel, real scalar islog,
    real scalar useconstant, real scalar nosplit, real scalar zstar,
    real scalar t0, real scalar t1, real rowvector bcoef, real scalar elast,
    real scalar zlo, real scalar zhi, real scalar zL, real scalar zH,
    real scalar bw)
{
    real rowvector out
    out = pbx_bias_core_ws(estimator, bmodel, islog, useconstant, nosplit,
        zstar, t0, t1, bcoef, elast, zlo - zstar, zhi - zstar, zL - zstar,
        zH - zstar, bw, 1)
    // pbx_bias_core_ws echoes columns 13-16 as OFFSETS from the cutoff
    // (lo,hi,L,H, its own native convention); this pre-2.2.0 signature's
    // callers expect the ABSOLUTE zlo/zhi/zL/zH they passed in, matching
    // pb_fobias_core()/polbunchbias_legacy's convention exactly -- patch
    // them back (xscale=1 here, so no other column needs adjustment: dL/dR
    // are shift-invariant by construction and already match).
    out[13] = zlo
    out[14] = zhi
    out[15] = zL
    out[16] = zH
    return(out)
}

// -------------------------------------------------------------------
// Estimator-2 concentrated objective.  Given a TRIAL delc (Chetty's
// excess-mass ratio, as it enters the bias formula -- not the point
// estimator's own delta), solve the BETA-ONLY (K+1)-square normal
// equations for the beta-bias with delc held fixed (q=1/(1+delc) is then
// just a number, not a dimension of the system being inverted -- this is
// what structurally eliminates Step 12's "delta row/column collapses"
// failure mode, rather than merely detecting it after the fact), and
// return the resulting concentrated sum-of-squared-residuals Q(delc):
// the standard "y'y - Gtu'beta_hat" identity for a linear GLS/OLS problem,
// where y'y is the (right-region fit + mass-constraint) sum of squares
// BEFORE any beta correction.  betapar is written with the solved
// beta-bias (Mata's usual write-through-argument convention); missing Q
// (and missing betapar) signal an ill-conditioned trial delc -- the
// multi-start search treats these exactly like a bad grid point, never a
// crash.
// -------------------------------------------------------------------
real scalar pbx_bias_e2_objQ(real scalar delc, real scalar nosplit,
    real matrix Glo, real matrix Ghi, real rowvector Rlo, real rowvector Rhi,
    real rowvector Rbar, real scalar trueMass, real rowvector betat,
    real rowvector h1coef, real scalar H, real scalar hi,
    real colvector dsc, real colvector betapar)
{
    real scalar q, K
    real matrix GtGb
    real colvector Gtub
    real rowvector ures, Jmb
    real scalar uM, yTy

    q = 1/(1+delc)
    K = cols(betat) - 1

    ures = pbx_polsub(h1coef, q :* betat)

    GtGb = Glo + q^2 :* Ghi
    Gtub = q :* pbx_polmoment(ures, H, hi, K)

    // Excluded-region counterfactual in the Chetty mass row.
    //   splitmass: h0 below z*, h1 = h0/(1+delc) above z*  -> Rlo + q*Rhi
    //   poolmass : whole region at h0 (no deflation)        -> Rlo + Rhi
    // Both then add the delc/(1+delc) int_{z*}^{zbar} h0 missing-mass term.
    if (nosplit) Jmb = Rlo + Rhi + (delc*q) :* Rbar
    else         Jmb = Rlo + q :* Rhi + (delc*q) :* Rbar
    uM  = trueMass - (Jmb * betat')

    GtGb = GtGb + Jmb' * Jmb
    Gtub = Gtub + Jmb' * uM

    if (pbx_illcond((dsc * dsc') :* GtGb)) {
        betapar = J(K+1, 1, .)
        return(.)
    }

    betapar = dsc :* (invsym((dsc * dsc') :* GtGb) * (dsc :* Gtub))
    // y'y is the RAW target sum of squares (ures on the right window,
    // uM for the mass row) -- the q that scales the DESIGN (P -> q*P) does
    // not also rescale the target itself.  Multiplying this term by q^2, as
    // an earlier version of this line did, breaks the y'y - Gtu'beta_hat
    // identity (which must be >=0 for any well-posed least squares fit)
    // and lets Q run away to large negative values as q shrinks.
    yTy = pbx_polint(pbx_polmul(ures, ures), H, hi) + uM^2
    return(yTy - (Gtub' * betapar))
}

// -------------------------------------------------------------------
// Estimator-2 fixed-point iteration for ONE starting value delc0.  This
// is EXACTLY the old single-start Gauss-Newton loop's per-iteration body
// (same (K+2)-square GtG/Gtu assembly, same Mhi/qff cross terms built
// from the FIXED true coefficients betat, same 0.7-damped step, same
// convergence tolerance) -- deliberately NOT reformulated, because that
// Jacobian (which holds betat fixed rather than re-differentiating at
// the running bias estimate) defines a specific, already-validated fixed
// point equation, not merely an approximation to some other objective.
// An earlier version of this file replaced the equation itself (profiling
// beta out and globally minimizing the resulting concentrated SSER via
// golden section) -- which is a well-posed problem in its own right, but
// a DIFFERENT one: validating it against test_pbx_est2.do's independent
// K=1 reference (pb_fobias_core) showed differences up to 0.10, i.e. it
// silently changed the estimand, not just its robustness.  This function
// restores the original equation; only pbx_bias_e2_solve() below is new,
// calling this from many starting points instead of one.
// -------------------------------------------------------------------
void pbx_bias_e2_fixedpoint(real scalar delc0, real scalar nosplit,
    real matrix Glo, real matrix Ghi, real matrix Mhi, real scalar qff,
    real rowvector Rlo, real rowvector Rhi, real rowvector Rbar,
    real scalar Sbar, real scalar Sright, real scalar trueMass,
    real rowvector betat, real rowvector h1coef, real scalar H,
    real scalar hi, real colvector dsc2, real scalar K, real scalar K1,
    real scalar np, real scalar delc, real colvector biaspar,
    real scalar converged)
{
    real scalar q, dqdD, dstep, nlsit
    real matrix GtG
    real colvector Gtu
    real rowvector ures, Jm
    real scalar uM

    delc = delc0
    for (nlsit = 1; nlsit <= 80; nlsit++) {
        q    = 1/(1+delc)
        dqdD = -1/((1+delc)^2)

        ures = pbx_polsub(h1coef, q :* betat)

        GtG = J(np, np, 0)
        Gtu = J(np, 1, 0)

        GtG[|1,1 \ K1,K1|] = Glo + q^2 :* Ghi
        GtG[|1,np \ K1,np|] = (q*dqdD) :* Mhi
        GtG[|np,1 \ np,K1|] = ((q*dqdD) :* Mhi)'
        GtG[np,np] = dqdD^2 * qff

        Gtu[|1 \ K1|] = q :* pbx_polmoment(ures, H, hi, K)
        Gtu[np] = dqdD * pbx_polint(pbx_polmul(betat, ures), H, hi)

        Jm = J(1, np, 0)
        if (nosplit) {
            // poolmass: excluded-region cf is Rlo + Rhi (undeflated), so
            // the only delc-dependent mass-row piece is (delc*q)*Rbar.
            Jm[|1,1 \ 1,K1|] = Rlo + Rhi + (delc*q) :* Rbar
            Jm[1,np]         = (q^2)*Sbar
            uM = trueMass - ((Rlo + Rhi + (delc*q) :* Rbar) * betat')
        }
        else {
            Jm[|1,1 \ 1,K1|] = Rlo + q :* Rhi + (delc*q) :* Rbar
            Jm[1,np]         = dqdD*Sright + (q^2)*Sbar
            uM = trueMass - ((Rlo + q :* Rhi + (delc*q) :* Rbar) * betat')
        }

        GtG = GtG + Jm' * Jm
        Gtu = Gtu + Jm' * uM

        // No per-iteration conditioning check here, matching the original
        // loop exactly: only the FINAL matrix is tested below (Step 12's
        // fix).  A transient bad-conditioning step mid-iteration that the
        // trajectory later recovers from must not abort early -- the old
        // loop never looked at intermediate steps either.
        biaspar = dsc2 :* (invsym((dsc2 * dsc2') :* GtG) * (dsc2 :* Gtu))

        dstep = 0.7 * biaspar[np]
        delc  = delc + dstep
        // clamp, don't abort -- matches the original loop exactly, since
        // a transient dip past this boundary can still recover on a later
        // iteration (confirmed by the window-width sweep this was
        // validated against: aborting here instead of clamping silently
        // turned several recoverable widths into all-missing).
        if (1 + delc <= 1e-6) delc = -1 + 1e-6
        if (abs(dstep) < 1e-11) break
    }

    // post-loop conditioning check on the matrix from the last iteration
    // actually solved (Step 12's fix): a collapsed delta row/column here
    // means this starting value's trajectory landed somewhere delta isn't
    // identified, and the fixed point is not to be trusted.
    if (pbx_illcond((dsc2 * dsc2') :* GtG)) {
        delc = .
        biaspar = J(np, 1, .)
        converged = 0
        return
    }
    converged = 1
}

// -------------------------------------------------------------------
// Estimator-2 multi-start delta-solve: runs pbx_bias_e2_fixedpoint() from
// several starting values instead of the old single start at delc=Delta,
// mirroring polbunch.ado's own prof_delta_solve for the point-estimate
// profile.  A single starting point is exactly the failure mode
// prof_delta_solve was built to fix (wrong-basin convergence); separately,
// the window-width sweep this was validated against showed the single
// path from Delta can wander through a badly-conditioned region as delc
// moves away from its start and settle on a wrong-but-smooth root as the
// window widens.  Other starting points often avoid that path entirely.
//
// When pbx_bias_e2_fixedpoint() converges from >1 starting value to
// DIFFERENT roots, pbx_bias_e2_objQ()'s concentrated sum-of-squared
// residuals (a well-defined merit function even though it is not itself
// what the Gauss-Newton Jacobian's fixed point solves for) breaks the
// tie by preferring the root that fits the data best.  Delta itself is
// always among the starting values, so whenever the fixed point is
// unique, this reduces to exactly the old single-start answer.
// -------------------------------------------------------------------
void pbx_bias_e2_solve(real scalar nosplit,
    real matrix Glo, real matrix Ghi, real matrix Mhi, real scalar qff,
    real rowvector Rlo, real rowvector Rhi, real rowvector Rbar,
    real scalar Sbar, real scalar Sright, real scalar trueMass,
    real rowvector betat, real rowvector h1coef, real scalar H,
    real scalar hi, real colvector dsc, real colvector dsc2,
    real scalar K, real scalar K1, real scalar np, real scalar Delta,
    real scalar bestdelc, real colvector bestbeta, real scalar bestQ,
    real scalar nbasin, real scalar gapQ)
{
    real colvector starts, roots, Qs, biaspar, bestbeta_k1
    real scalar ns, i, delc, converged, j, dupe, Q, b2Q

    starts = (Delta \ 0 \ -0.9 \ -0.5 \ 0.5 \ 1 \ 2 \ 5 \ 10 \ -0.99)
    ns = rows(starts)
    roots = J(0, 1, .)
    Qs    = J(0, 1, .)

    for (i = 1; i <= ns; i++) {
        pbx_bias_e2_fixedpoint(starts[i], nosplit, Glo, Ghi, Mhi, qff, Rlo, Rhi, Rbar,
            Sbar, Sright, trueMass, betat, h1coef, H, hi, dsc2, K, K1, np,
            delc, biaspar, converged)
        if (!converged | missing(delc)) continue

        dupe = 0
        for (j = 1; j <= rows(roots); j++) {
            if (abs(roots[j] - delc) < 1e-6) {
                dupe = 1
                break
            }
        }
        if (dupe) continue

        Q = pbx_bias_e2_objQ(delc, nosplit, Glo, Ghi, Rlo, Rhi, Rbar, trueMass, betat,
                h1coef, H, hi, dsc, bestbeta_k1)
        roots = roots \ delc
        Qs    = Qs \ Q
    }

    bestdelc = .
    bestQ    = .
    bestbeta = J(K1, 1, .)
    nbasin   = rows(roots)
    gapQ     = .
    if (rows(roots) == 0) return

    for (i = 1; i <= rows(roots); i++) {
        if (missing(Qs[i])) continue
        if (missing(bestQ) | Qs[i] < bestQ) {
            bestQ    = Qs[i]
            bestdelc = roots[i]
        }
    }
    if (missing(bestdelc)) {
        nbasin = 0
        return
    }

    if (nbasin >= 2) {
        b2Q = .
        for (i = 1; i <= rows(roots); i++) {
            if (roots[i] == bestdelc | missing(Qs[i])) continue
            if (missing(b2Q) | Qs[i] < b2Q) b2Q = Qs[i]
        }
        if (bestQ > 0 & !missing(b2Q)) gapQ = b2Q/bestQ
    }

    // capture the return value: an unassigned Mata function-call statement
    // auto-echoes its result to the results window, which otherwise leaked
    // this SSR value (twice, whenever pbx_bias_core_ws retries with the
    // same inputs after a missing first attempt) as a stray bare number.
    bestQ = pbx_bias_e2_objQ(bestdelc, nosplit, Glo, Ghi, Rlo, Rhi, Rbar, trueMass, betat,
        h1coef, H, hi, dsc, bestbeta)
}

// ===================================================================
real rowvector pbx_bias_core_ws(
    real scalar estimator,
    real scalar bmodel,
    real scalar islog,
    real scalar useconstant,
    real scalar nosplit,
    real scalar zstarw,
    real scalar t0,
    real scalar t1,
    real rowvector bcoef,
    real scalar elast,
    real scalar lo,
    real scalar hi,
    real scalar L,
    real scalar H,
    real scalar bw,
    real scalar xscale
)
{
    real scalar tau, Ltau, x, rho, rhow, Delta, r, B, K, a0, m0
    real scalar dL, dR, lambda_rep, pref
    real scalar bias_h, bias_B, bias_resp, bias_shift, bias_e
    real scalar bias_slope, bias_lambda
    real scalar atilde, mtilde, Btilde, rtilde, rr, guess
    real scalar edge_ovh, overhang
    real scalar zstar_out, lo_out, hi_out, L_out, H_out
    real matrix M, GtG, Glo, Ghi
    real colvector Gtu, Rvec, biaspar, Mhi, dsc, dsc2, bbias
    real rowvector ucoef, bhat, h1coef, betat, ures, Rlo, Rhi, Rbar, Jm
    real scalar hminus, hplus, Hstar, Bsaez, sleft, sright
    real scalar sright0, hright0, m_saez, a_saez
    real scalar Asaez, qsaez, disc, xhat, dlogzhat
    real scalar K1, np, delc, q, dqdD, dstep, nlsit, qff, ddel
    real scalar uM, Sbar, Sright, lo_right, trueRightMass, trueMass, hlo, llo
    real scalar bestQ2, nbasin2, gapQ2
    real rowvector biasbeta

    tau   = (1-t0)/(1-t1)
    Ltau  = ln(tau)
    x     = tau^elast
    rho   = ln(x)
    // rho is a shift in the RUNNING VARIABLE's own units (log mode's
    // analogue of the level case's z*.Delta level shift) -- w-space
    // formulas (pbx_polshift arguments, additive offsets to w-space
    // window bounds) need it in w-units, exactly like zstar -> zstarw;
    // rho itself (columns 11 / Delta's sibling) is reported unconverted,
    // same as x/Delta/tau, since it is a pure function of t0/t1/elast.
    rhow  = rho/xscale
    Delta = x - 1
    if (bw <= 0) bw = 1
    if (xscale <= 0) xscale = 1

    K  = cols(bcoef) - 1
    a0 = bcoef[1]
    m0 = 0
    if (K >= 1) m0 = bcoef[2]

    // z-unit echoes for display only (see header note); lo,hi,L,H,zstarw
    // themselves are used, unmodified, in every formula below
    zstar_out = zstarw*xscale
    lo_out    = lo*xscale
    hi_out    = hi*xscale
    L_out     = L*xscale
    H_out     = H*xscale

    if (islog) r = rhow
    else       r = zstarw*Delta

    B = pbx_polint(bcoef, 0, r)/bw

    // dL/dR are shift-invariant (their z/zL/zlo coefficients sum to 0), so
    // the cutoff's ABSOLUTE position can be set to 0 here (w-space is
    // already centred there) without needing zstarw at all; convert the
    // w-space result to a z-unit echo at the end, same as lo_out etc.
    dL = (0 - L + (L-lo)/2) * xscale
    dR = (H - 0 + (hi-H)/2) * xscale

    lambda_rep = pbx_lambda_report(bcoef, zstarw, islog, xscale)

    bias_h = bias_B = bias_resp = bias_shift = bias_e = .
    bias_slope = bias_lambda = .
    biasbeta = J(1, K + 1, .)

    if (estimator == 1) {
        // ---- fitting-stage bias -----------------------------------------
        // plim(betahat) - beta = (int_W P P')^{-1} int_right P (h1_true - h0)
        // W = [lo,L] u [H,hi] ; contamination is nonzero only on the right.
        // diagonal-preconditioned normal equations (exact change of basis;
        // keeps a wide centred window from tripping pbx_illcond spuriously)
        dsc = pbx_monoscale(lo, hi, L, H, K)
        M   = (dsc * dsc') :* (pbx_polgram(lo, L, K) + pbx_polgram(H, hi, K))
        if (pbx_illcond(M)) {
            return((estimator, bmodel, islog, zstar_out, t0, t1, tau, lambda_rep,
                elast, x, rho, Delta, lo_out, hi_out, L_out, H_out, dL, dR, B,
                J(1, 7, .), J(1, K+1, .)))
        }

        // NOTE (log mode, 2026-09-06): investigated and largely CLEARED.
        // The log branch here matches polbunch (h1 = h0(w+rhow), the
        // additive iso-elastic shift) and is byte-exact vs a brute-force
        // grid.  Oracle plim checks with the true h0 reproduce the plim
        // naive elasticity to ~1e-3 in BOTH level and log, both `constant'
        // and `exact' inversion, for estimators 1 and 2 -- i.e. the
        // log-mode bias formula is CORRECT.  Remaining minor item: the
        // estimator-2 Newton `iterate' shape-solve is a little less robust
        // to polynomial OVER-fitting in log mode than in level (a degree-5
        // fit of a degree-3 counterfactual leaves ~0.04 in log vs ~0.001
        // in level); at matched degree log mode lands on the truth.
        if (islog) ucoef = pbx_polsub(pbx_polshift(bcoef, rhow), bcoef)
        else       ucoef = pbx_polsub(pbx_polaffine(bcoef, x, r), bcoef)

        Gtu     = pbx_polmoment(ucoef, H, hi, K)
        biaspar = dsc :* (invsym(M) * (dsc :* Gtu))

        // plausibility guard.  The fitting-stage bias is a first-order
        // (small-contamination) object and everything downstream -- the
        // linearised relative-slope bias, the response inversion -- assumes
        // ||biaspar|| is small next to ||bcoef||.  On a wide centred window at
        // a high polynomial order the naive one-sided monomial design can be
        // collinear enough to slip past pbx_illcond() yet still return a
        // coefficient bias that dwarfs the counterfactual it corrects (a bias
        // in h0(zstar) larger than h0(zstar) itself).  The reported numbers
        // are then a numerical artifact, not an economic bias -- bail exactly
        // as the ill-conditioned branch above does.
        if (abs(a0) <= 1e-14 | abs(biaspar[1]) > abs(a0)) {
            return((estimator, bmodel, islog, zstar_out, t0, t1, tau, lambda_rep,
                elast, x, rho, Delta, lo_out, hi_out, L_out, H_out, dL, dR, B,
                J(1, 7, .), J(1, K+1, .)))
        }

        biasbeta = biaspar[|1 \ K+1|]'
        bias_h = biaspar[1]
        if (K >= 1) bias_slope = biaspar[2]

        // reported relative-slope bias -- linearised, exactly as
        // the local-linear polbunchbias.  Pure diagnostic: feeds nothing downstream.
        // (An exact ratio-difference version lands with the .ado rewrite.)
        if (K >= 1 & abs(a0) > 1e-14) {
            pref = zstarw
            if (islog) pref = 1/xscale
            bias_lambda = pref * (bias_slope/a0 - (m0/a0)*(bias_h/a0))
        }

        // ---- bunching-mass bias ----------------------------------------
        // fit-error piece: -(1/bw) int_E (h0hat - h0) = -(1/bw) int_E P'dbeta
        Rvec   = pbx_polmoment((1), L, H, K)
        bias_B = -(Rvec' * biaspar) / bw

        // overhang / compression term (upper edge only; needs L<0<H)
        if (L < 0 & H > 0) {
            if (islog) edge_ovh = H + rhow
            else       edge_ovh = x*H + r
            overhang = pbx_polint(bcoef, H, edge_ovh)/bw - B
            bias_B   = bias_B + overhang
        }

        // ---- response inversion --------------------------------------
        bhat   = bcoef + biaspar'
        atilde = bhat[1]
        Btilde = B + bias_B

        rtilde = .
        if (abs(atilde) > 1e-14) rtilde = Btilde*bw/atilde   // constant-density
        if (!useconstant) {
            // exact: invert the fitted-density integral; if that fails
            // (noisy estimate with no usable root) keep the constant value
            rr = pbx_polroot(bhat, Btilde*bw, rtilde)
            if (!missing(rr)) rtilde = rr
        }

        bias_resp = rtilde - r
        if (islog) {
            bias_shift = .
            bias_e     = bias_resp*xscale/Ltau
        }
        else {
            bias_shift = rtilde/zstarw - Delta
            bias_e     = ln(1 + rtilde/zstarw)/Ltau - elast
        }
    }

    if (estimator == 2) {
        // Chetty-style estimator: fits h0hat = P'beta on the left, and
        // h1hat = q*(P'beta) with q = 1/(1+delta) on the right, plus the
        // Chetty mass-row constraint
        //   Hstar = int_E h0hat_left + q*int_E h0hat_right
        //           + delta/(1+delta) * int_{zstar}^{zbar} h0hat.
        // The missing-mass term carries delta/(1+delta) (= delta*int h1),
        // NOT delta*int h0, because h1 = h0/(1+delta); this matches
        // bmodel_row23()/bmodel_row() in polbunch.ado.  NLS in (beta, delta), solved by a
        // Gauss-Newton fixed point on delta (delc) linearised about the
        // truth.  The (K+2)-square system generalises the local-linear
        // 3x3.  Structure is unit-agnostic; only the right-side truth
        // (h1coef) and trueMass differ between levels and logs.
        betat = bcoef
        K1 = K + 1
        np = K + 2

        if (islog) h1coef = pbx_polshift(betat, rhow)
        else       h1coef = pbx_polaffine(betat, x, r)

        // delta-independent mass-row window integrals (all divided by bw)
        Rlo = J(1, K1, 0)
        if (L < 0) {
            hlo = H
            if (H > 0) hlo = 0
            Rlo = pbx_polmoment((1), L, hlo, K)' / bw
        }
        Rhi = J(1, K1, 0)
        if (H > 0) {
            llo = L
            if (L < 0) llo = 0
            Rhi = pbx_polmoment((1), llo, H, K)' / bw
        }
        Rbar   = pbx_polmoment((1), 0, hi, K)' / bw
        Sbar   = (Rbar * betat')
        Sright = (Rhi  * betat')

        // true observed mass in the excluded window [L,H]
        if (islog) {
            trueMass = pbx_polint(betat, L, H + rhow)/bw
        }
        else {
            lo_right = 0
            if (L > 0) lo_right = L
            trueRightMass = 0
            if (H > lo_right) {
                trueRightMass = ///
                    pbx_polint(betat, x*lo_right + r, x*H + r)/bw
            }
            trueMass = (Rlo * betat') + trueRightMass
            if (L < 0 & H > 0) trueMass = trueMass + B
        }

        // delta-independent monomial integrals over the fitting windows
        Glo = pbx_polgram(lo, L, K)              // left P P'
        Ghi = pbx_polgram(H,  hi, K)             // right P P'
        Mhi = pbx_polmoment(betat, H, hi, K)     // right P (P'betat)
        qff = (betat * Ghi * betat')             // right (P'betat)^2

        // diagonal preconditioner: sc^{-k} on the K+1 beta coefficients
        // (monoscale, normalises the s-power spread).  The delta index
        // needs its OWN factor: the delta row/column of the (K+2)-square
        // Gauss-Newton system scales with the overall HEIGHT of h0 -- Mhi
        // is linear in betat, qff = betat*Ghi*betat' is quadratic -- while
        // the beta block is h0-scale-invariant, so a tall h0 (counts ~1e5,
        // or an h0poly() given in count units) leaves GtG[np,np] ~ height^2
        // dominating the O(1) beta block and trips pbx_illcond even when
        // the system is perfectly solvable.  Scaling the delta index by
        // 1/sqrt(qff) sends GtG[np,np] -> O(dqdD^2) and the off-diagonal
        // -> O(1).  Still an exact change of basis (the D's around
        // invsym(D GtG D) cancel), so it moves only floating-point paths
        // and which cases clear the conditioning guard -- never a
        // well-posed result.
        dsc  = pbx_monoscale(lo, hi, L, H, K)
        ddel = 1
        if (qff > 1e-300) ddel = 1/sqrt(qff)
        dsc2 = (dsc \ ddel)
        if (pbx_illcond((dsc * dsc') :* (Glo + Ghi))) {
            return((estimator, bmodel, islog, zstar_out, t0, t1, tau, lambda_rep,
                elast, x, rho, Delta, lo_out, hi_out, L_out, H_out, dL, dR, B,
                J(1, 7, .), J(1, K+1, .)))
        }

        // Multi-start delta-solve (STEP 14): runs the SAME (K+2)-square
        // Gauss-Newton fixed point as before from several starting values
        // of delc instead of the single start at delc=Delta -- the same
        // wrong-basin risk prof_delta_solve was built to fix for
        // polbunch's own point estimator.  (An earlier version of this
        // step reformulated the equation itself as a globally-minimized
        // concentrated SSR; that is a well-posed problem but a DIFFERENT
        // one from what this Jacobian's fixed point solves for --
        // validating against test_pbx_est2.do's independent K=1 reference
        // showed differences up to 0.10, i.e. it silently changed the
        // estimand.  See pbx_bias_e2_fixedpoint()'s header.)  Delta is
        // always among the starting values, so a unique fixed point
        // reproduces the old single-start answer exactly.  bestQ2/
        // nbasin2/gapQ2 are computed (weak-ID / competing-basin
        // diagnostics, as in prof_delta_solve) but not yet surfaced to
        // r()/e() -- a deliberate first-pass scope limit, not an
        // oversight.
        pbx_bias_e2_solve(nosplit, Glo, Ghi, Mhi, qff, Rlo, Rhi, Rbar, Sbar, Sright,
            trueMass, betat, h1coef, H, hi, dsc, dsc2, K, K1, np, Delta,
            delc, bbias, bestQ2, nbasin2, gapQ2)
        if (missing(delc)) {
            return((estimator, bmodel, islog, zstar_out, t0, t1, tau, lambda_rep,
                elast, x, rho, Delta, lo_out, hi_out, L_out, H_out, dL, dR, B,
                J(1, 7, .), J(1, K+1, .)))
        }
        biaspar = bbias

        // plausibility guard (see the estimator-1 branch): a coefficient bias
        // that dwarfs the counterfactual it corrects is a numerical artifact
        // of a near-collinear design, not an economic bias.
        if (abs(a0) <= 1e-14 | abs(biaspar[1]) > abs(a0)) {
            return((estimator, bmodel, islog, zstar_out, t0, t1, tau, lambda_rep,
                elast, x, rho, Delta, lo_out, hi_out, L_out, H_out, dL, dR, B,
                J(1, 7, .), J(1, K+1, .)))
        }

        biasbeta = biaspar[|1 \ K1|]'
        bias_h = biaspar[1]
        if (K >= 1) bias_slope = biaspar[2]
        if (K >= 1 & abs(a0) > 1e-14) {
            pref = zstarw
            if (islog) pref = 1/xscale
            bias_lambda = pref * (bias_slope/a0 - (m0/a0)*(bias_h/a0))
        }

        // The bmodel(1) command option was removed, so `bmodel' is always 0
        // from polbunchbias itself and the reduced-form Bhat branch below is
        // what estimator 2 always uses (governed by poolmass / splitmass).
        // The bmodel==1 branch is retained for direct Mata callers / the
        // regression fixtures only.
        if (bmodel == 1 & !islog) {
            // model-implied response: the NLS pseudo-true delta itself
            bias_shift = delc - Delta
            bias_resp  = zstarw*bias_shift
            bias_B     = .
            bias_e     = ln(1 + Delta + bias_shift)/Ltau - elast
        }
        else {
            // reduced-form Bhat = Hstar - int_E h0hat
            Rvec   = pbx_polmoment((1), L, H, K)
            bias_B = -(Rvec' * biaspar[|1 \ K1|]) / bw
            bhat   = betat + biaspar[|1 \ K1|]'

            if (L < 0 & H > 0) {
                if (islog) edge_ovh = H + rhow
                else       edge_ovh = x*H + r
                overhang = pbx_polint(betat, H, edge_ovh)/bw - B
                bias_B = bias_B + overhang
                // splitmass: polbunch forms Bhat as
                //   Hstar - int_{zL}^{z*} h0hat
                //         - int_{z*}^{zH} h0hat/(1+deltahat)
                // i.e. it deltahat-deflates the FITTED h0hat above the kink
                // (the Chetty h1 = h0/(1+delta) restriction), NOT an
                // iso-elastic (1-1/x) stretch of the TRUE h0.  Following
                // that procedure, the extra term over poolmass is
                //   (deltahat/(1+deltahat)) * int_{z*}^{zH} h0hat / bw
                // with deltahat = delc (the NLS pseudo-true delta solved
                // for above) and h0hat = bhat.  See test_pbx_est2.do B1.
                if (!nosplit)
                    bias_B = bias_B ///
                        + (delc/(1 + delc))*pbx_polint(bhat, 0, H)/bw
            }

            atilde = bhat[1]
            Btilde = B + bias_B

            rtilde = .
            if (abs(atilde) > 1e-14) rtilde = Btilde*bw/atilde
            if (!useconstant) {
                rr = pbx_polroot(bhat, Btilde*bw, rtilde)
                if (!missing(rr)) rtilde = rr
            }

            bias_resp = rtilde - r
            if (islog) {
                bias_shift = .
                bias_e     = bias_resp*xscale/Ltau
            }
            else {
                bias_shift = rtilde/zstarw - Delta
                bias_e     = ln(1 + rtilde/zstarw)/Ltau - elast
            }
        }
    }

    if (estimator == 3 | estimator == 0) {
        // Estimators 3 and 0 have NO fitting-stage bias at any degree K:
        // (beta) are the truth (est 3 imposes the correct shift-and-stretch
        // restriction; est 0 fits h0 from uncontaminated left data only).
        // Their only biases are the two downstream researcher choices --
        // nosplit (overhang) and useconstant (density shortcut).
        bias_h = 0
        bias_slope = 0
        bias_lambda = 0
        biasbeta = J(1, K + 1, 0)

        bias_B = 0
        if (nosplit & L < 0 & H > 0) {
            if (islog) edge_ovh = H + rhow
            else       edge_ovh = x*H + r
            bias_B = pbx_polint(bcoef, H, edge_ovh)/bw - B
        }

        Btilde = B + bias_B
        rtilde = .
        if (abs(a0) > 1e-14) rtilde = Btilde*bw/a0        // constant-density
        if (!useconstant) {
            // invert the exact integral of the (uncontaminated) h0
            rr = pbx_polroot(bcoef, Btilde*bw, rtilde)
            if (!missing(rr)) rtilde = rr
        }

        if (bias_B == 0 & !useconstant) {
            bias_resp  = 0
            bias_shift = 0
            bias_e     = 0
        }
        else {
            bias_resp = rtilde - r
            if (islog) {
                bias_shift = .
                bias_e     = bias_resp*xscale/Ltau
            }
            else {
                bias_shift = rtilde/zstarw - Delta
                bias_e     = ln(1 + rtilde/zstarw)/Ltau - elast
            }
        }
    }

    if (estimator == 4) {
        // Saez three-region trapezoid.  Intrinsically a TWO-POINT method:
        // it recovers only a line, so its bias against a curved h0 comes
        // from evaluating the exact (curved) h0 over the reference regions.
        // The mass->response transform equation is unchanged -- Saez only
        // ever knows the two reference densities hminus, hplus.
        if (islog) h1coef = pbx_polshift(bcoef, rhow)
        else       h1coef = pbx_polaffine(bcoef, x, r)

        // reference densities = exact region MEANS (not midpoint values)
        hminus = pbx_polint(bcoef,  lo, L )/(L  - lo)
        hplus  = pbx_polint(h1coef, H,  hi)/(hi - H )

        // true observed mass in the excluded window [L,H]:
        //   h0 below the cutoff, point mass B at it, h1_true above it
        Hstar = pbx_polint(bcoef, L, 0)/bw + B + pbx_polint(h1coef, 0, H)/bw
        Bsaez = Hstar - (-L/bw)*hminus - (H/bw)*hplus

        // Saez's implicit two-point counterfactual line through the exact
        // h0 points (sleft, hminus) and (sright0, hright0)
        sleft  = (lo + L)/2
        sright = (H  + hi)/2
        if (islog) {
            sright0 = sright + rhow
            hright0 = hplus
        }
        else {
            sright0 = (x - 1)*zstarw + x*sright
            hright0 = hplus/x
        }
        m_saez = (hright0 - hminus)/(sright0 - sleft)
        a_saez = hminus - m_saez*sleft

        bias_h     = a_saez - a0
        bias_slope = m_saez - m0
        biasbeta = J(1, K + 1, 0)
        biasbeta[1] = bias_h
        if (K >= 1) biasbeta[2] = bias_slope
        pref = zstarw
        if (islog) pref = 1/xscale
        bias_lambda = pref*(m_saez/a_saez - m0/a0)

        if (islog) {
            if (hminus > 0 & hplus > 0) {
                if (useconstant) dlogzhat = Bsaez*bw/hminus
                else             dlogzhat = 2*Bsaez*bw/(hminus + hplus)
                xhat       = exp(dlogzhat)
                bias_B     = Bsaez - B
                bias_resp  = dlogzhat - rhow
                bias_shift = .
                bias_e     = bias_resp*xscale/Ltau
            }
        }
        else {
            if (hminus > 0) {
                Asaez = 2*Bsaez*bw/zstarw
                if (useconstant) {
                    xhat       = 1 + Asaez/(2*hminus)
                    bias_B     = Bsaez - B
                    bias_resp  = zstarw*(xhat - x)
                    bias_shift = xhat - x
                    bias_e     = ln(xhat)/Ltau - elast
                }
                else {
                    qsaez = hplus - hminus - Asaez
                    disc  = qsaez^2 + 4*hminus*hplus
                    if (disc >= 0) {
                        xhat = (-qsaez + sqrt(disc))/(2*hminus)
                        if (xhat > 0) {
                            bias_B     = Bsaez - B
                            bias_resp  = zstarw*(xhat - x)
                            bias_shift = xhat - x
                            bias_e     = ln(xhat)/Ltau - elast
                        }
                    }
                }
            }
        }
    }

    // ---- convert the w-space-only quantities back to z-units --------
    // bias_h, bias_B, bias_shift, bias_e, bias_lambda are ALREADY correct
    // (coordinate-invariant given zstarw/bw used consistently throughout
    // -- see the header note).  bias_resp is a LEVEL shift (needs *xscale)
    // and bias_slope is a derivative (needs /xscale); biasbeta is the
    // per-coefficient vector, so pbx_polscale(.,1/xscale) applies the
    // matching 1/xscale^k to each entry k in one step.  Multiplying/
    // dividing a missing value by a finite xscale stays missing, so this
    // is safe to apply unconditionally, including on the all-missing
    // early-return branches above.
    if (!missing(bias_resp))  bias_resp  = bias_resp * xscale
    if (!missing(bias_slope)) bias_slope = bias_slope / xscale
    biasbeta = pbx_polscale(biasbeta, 1/xscale)

    // columns 1..26 (the local-linear layout); then the per-coefficient bias
    // vector biasbeta (length K+1), so the caller can recover the full
    // fitted polynomial bhat = bcoef + biasbeta.
    return((estimator, bmodel, islog, zstar_out, t0, t1, tau, lambda_rep, elast,
        x, rho, Delta, lo_out, hi_out, L_out, H_out, dL, dR, B, bias_h, bias_B, bias_resp,
        bias_shift, bias_e, bias_slope, bias_lambda, biasbeta))
}

// -------------------------------------------------------------------
// pbx_bias_scres_ws(): the self-consistency residual F(x) - xhat, where
// x = (h0 shape coeffs 2..K1 in the dsc-preconditioned basis, elasticity).
// A root of this is a counterfactual (h0 shape, elasticity) that, run
// through the bias map, reproduces the observed fit -- i.e. the fixed
// point the `iterate' option targets.  Returns all-missing on a
// singular design / non-invertible response.
// -------------------------------------------------------------------
real rowvector pbx_bias_scres_ws(real rowvector x, real rowvector tgt,
    real scalar estimator, real scalar bmodel, real scalar islog,
    real scalar useconstant, real scalar nosplit, real scalar zstarw,
    real scalar t0, real scalar t1, real scalar lo, real scalar hi,
    real scalar L, real scalar H, real scalar bw, real scalar xscale,
    real rowvector dsc)
{
    real scalar n, K1, e, lead
    real rowvector s, bp, b, res, bbeta, bpred, out

    n  = cols(x)
    K1 = n
    s  = x[|1 \ n - 1|]
    e  = x[n]
    bp = (1, s)
    b  = bp :* dsc
    res = pbx_bias_core_ws(estimator, bmodel, islog, useconstant, nosplit,
              zstarw, t0, t1, b, e, lo, hi, L, H, bw, xscale)
    out = J(1, n, .)
    if (missing(res[20]) | missing(res[24])) return(out)
    // biasbeta comes back /xscale^k unit-converted for reporting; undo
    // that so this residual stays in the native w-space basis dsc/bp use.
    bbeta = pbx_polscale(res[|27 \ 26 + K1|], xscale)
    bpred = bp + (bbeta :/ dsc)
    lead  = bpred[1]
    if (abs(lead) < 1e-14) return(out)
    bpred = bpred :/ lead
    out = (bpred[|2 \ K1|] :- tgt[|1 \ K1 - 1|], (e + res[24]) - tgt[K1])
    return(out)
}

// -------------------------------------------------------------------
// pbx_bias_newton_ws(): single-start damped Newton (finite-difference
// Jacobian + backtracking line search) for the pbx_bias_scres_ws root,
// started from (fitted h0 shape, naive elasticity).  Estimator 2's
// self-consistency fixed point is unique (verified by MC and on the
// Chetty fig-20a data) but the damped fixed-point substitution in
// pbx_bias_solve_ws is not a contraction there and runs away; Newton
// reaches it directly.  Returns 1 and fills bout/eout/niter on
// convergence to an admissible root, 0 otherwise (caller then falls
// back to the damped loop).
// -------------------------------------------------------------------
real scalar pbx_bias_newton_ws(
    real scalar estimator, real scalar bmodel, real scalar islog,
    real scalar useconstant, real scalar nosplit, real scalar zstarw,
    real scalar t0, real scalar t1, real rowvector bhat, real scalar ehat,
    real scalar lo, real scalar hi, real scalar L, real scalar H,
    real scalar bw, real scalar xscale, real rowvector dsc,
    real scalar allowneg, real scalar checkh0, real scalar EBOUND,
    real rowvector bout, real scalar eout, real scalar niter)
{
    real scalar n, K1, it, i, hstep, rn, alpha, ls, okj, done, mv
    real rowvector bhp, shp, tgt, x, Fx, xp, xm, Fp, Fm, xt, Rt, bp
    real matrix Jm
    real colvector dxv

    K1 = cols(bhat)
    n  = K1
    niter = 0
    bhp = bhat :/ dsc
    if (abs(bhp[1]) < 1e-14) return(0)
    bhp = bhp :/ bhp[1]
    shp = bhp[|2 \ K1|]
    tgt = (shp, ehat)
    x   = (shp, ehat)                 // single start

    for (it = 1; it <= 50; it++) {
        niter = it
        Fx = pbx_bias_scres_ws(x, tgt, estimator, bmodel, islog, useconstant,
                 nosplit, zstarw, t0, t1, lo, hi, L, H, bw, xscale, dsc)
        if (hasmissing(Fx)) return(0)
        rn = sqrt(Fx * Fx')
        if (rn < 1e-9) break
        Jm  = J(n, n, 0)
        okj = 1
        for (i = 1; i <= n; i++) {
            hstep = 1e-6 * max((abs(x[i]), 1e-3))
            xp = x ; xp[i] = xp[i] + hstep
            xm = x ; xm[i] = xm[i] - hstep
            Fp = pbx_bias_scres_ws(xp, tgt, estimator, bmodel, islog,
                     useconstant, nosplit, zstarw, t0, t1, lo, hi, L, H,
                     bw, xscale, dsc)
            Fm = pbx_bias_scres_ws(xm, tgt, estimator, bmodel, islog,
                     useconstant, nosplit, zstarw, t0, t1, lo, hi, L, H,
                     bw, xscale, dsc)
            if (hasmissing(Fp) | hasmissing(Fm)) okj = 0
            if (okj) Jm[, i] = ((Fp - Fm) / (2 * hstep))'
        }
        if (!okj) return(0)
        dxv = lusolve(Jm, Fx')
        if (hasmissing(dxv)) return(0)
        alpha = 1
        xt    = x
        done  = 0
        for (ls = 1; ls <= 25; ls++) {
            if (!done) {
                xt = x - alpha * dxv'
                Rt = pbx_bias_scres_ws(xt, tgt, estimator, bmodel, islog,
                         useconstant, nosplit, zstarw, t0, t1, lo, hi, L, H,
                         bw, xscale, dsc)
                if (!hasmissing(Rt) & sqrt(Rt * Rt') < rn) done = 1
                else alpha = alpha / 2
            }
        }
        if (!done) return(0)
        x = xt
    }

    Fx = pbx_bias_scres_ws(x, tgt, estimator, bmodel, islog, useconstant,
             nosplit, zstarw, t0, t1, lo, hi, L, H, bw, xscale, dsc)
    if (hasmissing(Fx) | sqrt(Fx * Fx') > 1e-6) return(0)

    eout = x[n]
    if (missing(eout) | abs(eout) > EBOUND) return(0)
    if (!allowneg & eout < 0) return(0)   // let the damped path apply the clamp
    bp = (1, x[|1 \ n - 1|])
    if (hasmissing(bp) | max(abs(bp)) > 1e6) return(0)
    bout = bp :* dsc
    if (checkh0) {
        mv = min((pbx_polminval(bout, lo, L), pbx_polminval(bout, H, hi)))
        if (mv < -1e-6) return(0)
    }
    return(1)
}

// -------------------------------------------------------------------
// pbx_bias_solve(): optional self-consistent iteration on
// (elasticity, h0 shape), then one final evaluation.  Writes results
// to <stub>_core (1x26), <stub>_bbeta (1xK1), <stub>_bcur (1xK1),
// <stub>_ecur, <stub>_iter, <stub>_conv.
// -------------------------------------------------------------------
void pbx_bias_solve_ws(
    real scalar estimator, real scalar bmodel, real scalar islog,
    real scalar useconstant, real scalar nosplit, real scalar zstarw,
    real scalar t0, real scalar t1, real rowvector bhat, real scalar ehat,
    real scalar lo, real scalar hi, real scalar L, real scalar H,
    real scalar bw, real scalar xscale, real scalar doiter, real scalar tol,
    real scalar maxiter, real scalar urlx, real scalar allowneg,
    real scalar checkh0, string scalar stub)
{
    real scalar K1, ecur, enew, iter, conv, dif, k
    real scalar diverged, EBOUND, newtok, enewt, nnewt
    real rowvector bcur, bnew, res, bbeta, bnewt
    real rowvector dsc, bcur_p, bnew_p, bbeta_p

    // |corrected elasticity| beyond this => the self-consistency fixed point
    // has run away (for a steep/curved level h0 the recentre-and-re-pin map
    // is not a contraction); the iterate is unusable and we fall back.
    EBOUND = 5

    K1 = cols(bhat)
    ecur = ehat
    bcur = bhat
    iter = 0
    conv = .
    diverged = 0

    // preconditioner: bhat is now a w-space (polbunch's own normalised
    // coordinate) coefficient vector, so it no longer spans many orders of
    // magnitude purely from polbunch's xscale the way a raw-s-coordinate
    // vector used to (that was the point of moving pbx_bias_core into
    // w-space).  Kept anyway as a cheap, exact change of basis -- a wide
    // w-space window (e.g. a long excluded region even in normalised
    // units) can still benefit -- and it makes this loop's own arithmetic
    // (differencing, re-pinning height, the >1e6 / dif<tol checks)
    // consistent with pbx_bias_core's internal Gram-solve preconditioning.
    // dsc[1] = 1 always, so height re-pinning is unaffected either way.
    dsc    = pbx_monoscale(lo, hi, L, H, K1 - 1)'
    bcur_p = bcur :/ dsc

    if (doiter) {
        conv   = 0
        newtok = 0
        // Estimators 1 and 2: solve the self-consistency fixed point with a
        // single-start Newton (finite-difference Jacobian + line search)
        // rather than the damped substitution.  For estimator 2 the damped
        // map is not a contraction and runs away (~100% at poly>=5); for
        // estimator 1 it converges only slowly.  The fixed point is unique
        // and sits at the truth (oracle plim checks; for est-2 splitmass
        // this needed the Chetty-deltahat harmonisation of bias_B first).
        // The damped loop stays as the fallback when Newton fails, and is
        // the only path for estimators 0/3/4.
        if (estimator == 1 | estimator == 2) {
            newtok = pbx_bias_newton_ws(estimator, bmodel, islog, useconstant,
                nosplit, zstarw, t0, t1, bhat, ehat, lo, hi, L, H, bw, xscale,
                dsc, allowneg, checkh0, EBOUND, bnewt, enewt, nnewt)
            if (newtok) {
                bcur   = bnewt
                bcur_p = bcur :/ dsc
                ecur   = enewt
                iter   = nnewt
                conv   = 1
            }
        }
        if (!newtok)
        for (k = 1; k <= maxiter; k++) {
            res    = pbx_bias_core_ws(estimator, bmodel, islog, useconstant,
                        nosplit, zstarw, t0, t1, bcur, ecur,
                        lo, hi, L, H, bw, xscale)
            bbeta  = res[|1, 27 \ 1, 26 + K1|]
            if (missing(res[20])) {
                iter = k
                break                        // singular design -- stop
            }
            // bbeta comes back from pbx_bias_core already unit-converted
            // (biasbeta /xscale^k) for external reporting; undo that here
            // so this loop's OWN arithmetic stays in the native w-space
            // basis bhat/bcur are already expressed in (pbx_polscale by
            // xscale is the exact inverse of the /xscale^k applied there).
            bbeta   = pbx_polscale(bbeta, xscale)
            bbeta_p = bbeta :/ dsc
            enew   = ehat - res[24]
            // economic floor: the compensated elasticity of taxable income
            // is non-negative (Slutsky) -- by default don't let the search
            // wander into a region no legitimate labour-supply model would
            // occupy (same logic as polbunch's delta>=0 default for the
            // convex-kink profile).  allowneg opts out.
            if (!allowneg & enew < 0) enew = 0
            bnew_p = (bhat :/ dsc) - bbeta_p
            if (abs(bnew_p[1]) > 1e-14) bnew_p = bnew_p :/ bnew_p[1]   // re-pin height
            // divergence guard: a blown-up coefficient vector or an
            // implausible corrected elasticity means this fixed point is
            // not a contraction here -- stop and fall back.  Optionally
            // (h0check, OFF by default) also reject a candidate that would
            // be NEGATIVE somewhere in the fitting region [lo,L] u [H,hi]
            // (same window pbx_illcond/pbx_monoscale treat as "the
            // relevant window"; the excluded gap [L,H] is pure structural
            // extrapolation, never checked).  Default OFF: an MC
            // validation against a known true bias found this check can
            // erase the entire benefit of promotion, not just filter bad
            // candidates -- h0 here is a low-degree polynomial
            // APPROXIMATION to a noisily estimated density, and such an
            // approximation can legitimately dip slightly negative in a
            // low-density region without the true h0 being invalid.
            if (missing(enew) | hasmissing(bnew_p) | max(abs(bnew_p)) > 1e6
                | abs(enew) > EBOUND
                | (checkh0 & min((pbx_polminval(bnew_p :* dsc, lo, L),
                    pbx_polminval(bnew_p :* dsc, H, hi))) < -1e-6)) {
                iter     = k
                diverged = 1
                break
            }
            dif    = max((abs(enew - ecur), rowmax(abs(bnew_p - bcur_p))))
            ecur   = ecur + urlx*(enew - ecur)
            bcur_p = bcur_p + urlx :* (bnew_p - bcur_p)
            bcur   = bcur_p :* dsc
            iter   = k
            if (dif < tol) {
                conv = 1
                break
            }
        }
        // a bounded run that merely exhausted maxiter without hitting `tol`
        // is KEPT (its last iterate is a partial refinement, not garbage);
        // only a runaway iterate is discarded.
        if (diverged) {
            bcur = bhat
            ecur = ehat
        }
    }

    res   = pbx_bias_core_ws(estimator, bmodel, islog, useconstant, nosplit,
                zstarw, t0, t1, bcur, ecur, lo, hi, L, H, bw, xscale)
    if (missing(res[20]) | missing(res[24]) |
        (doiter & abs(ehat - res[24]) > EBOUND)) {
        // the current point gives a missing or runaway bias -- fall back to
        // the raw (fitted) inputs so iterate never degrades the result
        if (doiter & (bcur != bhat | ecur != ehat)) diverged = 1
        bcur = bhat
        ecur = ehat
        conv = 0
        res  = pbx_bias_core_ws(estimator, bmodel, islog, useconstant, nosplit,
                    zstarw, t0, t1, bcur, ecur, lo, hi, L, H, bw, xscale)
    }
    bbeta = res[|1, 27 \ 1, 26 + K1|]

    // bcur is the accepted SELF-CONSISTENT counterfactual, stored for the
    // caller as r(corrected_beta); external convention is z-units (matches
    // bcoef's own convention, same as bbeta/biasbeta above), so convert out
    // of the native w-space basis this loop worked in, same as pbx_bias_core
    // does for biasbeta.
    st_matrix(stub + "_core",  res[|1, 1 \ 1, 26|])
    st_matrix(stub + "_bbeta", bbeta)
    st_matrix(stub + "_bcur",  pbx_polscale(bcur, 1/xscale))
    st_numscalar(stub + "_ecur", ecur)
    st_numscalar(stub + "_iter", iter)
    st_numscalar(stub + "_conv", conv)
    st_numscalar(stub + "_diverged", diverged)
}

end
