capture program drop polbunchbias
program define polbunchbias, rclass
    version 16.0

    /*
        iterate(): when to use it vs. not.

        The bias formula answers "if the true (elasticity, lambda) were X,
        how biased would this estimator be?" That's a different question
        depending on where X comes from:

        - Standalone mode with elasticity()/relslope() supplied from a KNOWN
          DGP (e.g. validating this formula against a Monte Carlo): X is
          already the truth. One evaluation of bias() at that point is the
          right answer; iterate() would search for a different fixed point
          than the one you already know is correct, so it should stay off
          (the default).

        - Post-estimation mode (letting polbunchbias pull elasticity/lambda
          from e()): both quantities come from the same biased estimator --
          lambda is computed from the fitted h0 polynomial, which is
          contaminated by the same bunching-induced misspecification that
          biases the elasticity estimate. Treating the raw estimates as a
          first guess at the truth and solving for the self-consistent
          fixed point theta_hat = theta_true + bias(theta_true) is exactly
          what iterate() does, which is why polbunch.ado's inline display
          turns it on by default.
    */
    syntax [, ESTimator(numlist max=1) ZSTAR(numlist max=1) ///
        T0(numlist max=1) T1(numlist max=1) ///
        RELSLOPE(numlist max=1) ELasticity(numlist max=1) ///
        ZLO(numlist max=1) ZHI(numlist max=1) ///
        ZL(numlist max=1) ZH(numlist max=1) ///
        BMODEL(integer 0) LOG ///
        BW(numlist max=1) ITERate TOLerance(real 1e-10) ///
        MAXITER(integer 100) UNDERRELAX(real 0.5) ///
        CONstant EXACT SPLITmass POOLmass ]

    /*
        Two toggle axes, each with a default that depends on estimator()
        (see bias.tex sec:threeaxes). Each is a pair of distinct flag
        names rather than an on/off flag, because Stata parses "no<flag>"
        as negation and that is indistinguishable from "unspecified" --
        which breaks an estimator-based default.

        Response inversion (step 3):
          constant : rhat = Bhat*bw / h0hat(z*)           (density shortcut)
          exact    : solve int_{z*}^{z*+rhat} h0hat = Bhat (integral)
          default  : constant for estimator(1)/(2), exact for estimator(0)/(3).

        Mass back-out (step 2):
          poolmass : Bhat = Hstar_obs - int_{zL}^{zH} h0hat -- pool the
                     fitted counterfactual over the whole window (Chetty
                     2011 eq. 16). Adds the overhang int_{z*}^{zH}(h1_true
                     - h0) to bias_B.
          splitmass: Bhat = Hstar_obs - int_{zL}^{z*} h0hat
                     - int_{z*}^{zH} h1hat -- split at the kink with the
                     estimator's own h1. No overhang.
          default  : poolmass for estimator(1)/(2), splitmass for
                     estimator(0)/(3). estimator(1) is forced to poolmass
                     (its h1 model IS h0).

        estimator(4) (Saez) reads only the response-inversion axis
        (constant = single side-band h-; exact = the two-point trapezoid);
        the mass axis does not apply.
    */
    if "`constant'" != "" & "`exact'" != "" {
        di as error "polbunchbias: specify at most one of constant / exact"
        exit 198
    }
    if "`poolmass'" != "" & "`splitmass'" != "" {
        di as error "polbunchbias: specify at most one of poolmass / splitmass"
        exit 198
    }

    // ----------------------------------------------------------------
    // Mode detection: count how many of the 9 primary required
    // options were supplied.  relslope is required too once in standalone
    // (explicit) mode, but it must NOT participate in this particular
    // count -- ncore==0 is what triggers e() mode, and "polbunchbias,
    // relslope(#)" alone (overriding just the slope, everything else read
    // from e()) has to keep landing there.
    // ----------------------------------------------------------------
    local ncore = ("`estimator'"!="") + ("`zstar'"!="") + ///
        ("`t0'"!="") + ("`t1'"!="") + ("`elasticity'"!="") + ///
        ("`zlo'"!="") + ("`zhi'"!="") + ("`zl'"!="") + ("`zh'"!="")

    if `ncore' == 0 {
        // ============================================================
        // e() MODE: extract all parameters from polbunch results
        // ============================================================
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

        /*
            Pull the two toggle axes from polbunch when it exposes them
            (e(nosplit), e(constant)); until it does, the estimator-based
            default below applies. An explicit flag on the command line
            still wins (resolved below).
        */
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
            di as error "(polbunch must be run with t0() and t1() to produce an elasticity estimate)"
            exit 111
        }
        local elasticity = `_elast'

        // lambda: compute from e(b) unless overridden by user via relslope()
        if "`relslope'" == "" {
            if `estimator' == 4 {
                // Saez: recover slope from the implicit two-point counterfactual
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
                local _den = `_sr0' - `_sl'
                if abs(`_den') < 1e-12 {
                    local lambda = 0
                }
                else {
                    local _m = (`_hr0' - `_hminus')/`_den'
                    local _a = `_hminus' - `_m'*`_sl'
                    local lambda = cond(`_a' > 0, `_m'*`zstar'/`_a', 0)
                }
            }
            else {
                // Polynomial: evaluate h0 and its derivative at zstar
                // using all polynomial terms in normalized coordinates.
                //
                // z_est = (z_orig - zmid)/xscale, so at z_orig = zstar:
                //   z_est = cutoff_est
                //   h0(zstar) = sum_k bk * cutoff_est^k
                //   dh0/dz_orig = (dh0/dz_est) / xscale
                //   lambda = (dh0/dz_est at cutoff_est)*zstar / (xscale*h0(cutoff_est))
                local _K      = e(polynomial)
                local _cest   = e(cutoff_est)
                local _xscale = e(xscale)
                local _h0c    = _b[h0:_cons]
                local _hval   = `_h0c'
                local _dhval  = 0
                local _zterm  "c.`_zname'"
                forvalues _k = 1/`_K' {
                    if `_k' > 1 local _zterm "`_zterm'#c.`_zname'"
                    capture local _bk = _b[h0:`_zterm']
                    if !_rc & !missing(`_bk') {
                        /*
                            cutoff_est is frequently negative (z gets
                            recentered around a midpoint, and the cutoff
                            often lands on the negative side). Substituting
                            a negative local straight into `_cest'^n lets
                            Stata parse the unary minus as binding looser
                            than "^", so -.12^0 becomes -(.12^0) = -1
                            instead of (-.12)^0 = 1 -- silently flipping the
                            sign of every term with an even exponent
                            (including 0). Parenthesize the base so the
                            whole negative value gets exponentiated, not
                            just its absolute value.
                        */
                        local _hval  = `_hval'  + `_bk' * (`_cest')^`_k'
                        local _dhval = `_dhval' + `_k' * `_bk' * (`_cest')^(`_k'-1)
                    }
                }
                if `_hval' > 0 {
                    local lambda = `_dhval' * `zstar' / (`_xscale' * `_hval')
                }
                else local lambda = 0
            }
        }
        else {
            local lambda = `relslope' * `zstar'
        }

        if "`bw'" == "" local bw = `_bw'
    }
    else if `ncore' == 9 & "`relslope'" != "" {
        // ============================================================
        // EXPLICIT MODE: all 9 primary options plus relslope() supplied
        // ============================================================
        local lambda = `relslope' * `zstar'
        if "`bw'"         == "" local bw = 1
        local islog_val = ("`log'" != "")
    }
    else {
        // Partial specification — error
        local missing_opts
        foreach v in estimator zstar t0 t1 elasticity zlo zhi zl zh relslope {
            if "``v''" == "" local missing_opts `missing_opts' `v'()
        }
        di as error "Specify all required options, or none to use polbunch results from e()"
        di as error "Missing: `missing_opts'"
        exit 198
    }

    if !inlist(`estimator', 0, 1, 2, 3, 4) {
        di as error "polbunchbias: only estimator(0), estimator(1), estimator(2), estimator(3), estimator(4) are supported"
        exit 198
    }

    /*
        Resolve the two toggle axes (see the option note above): explicit
        flag wins; else the e()-mode value if polbunch exposed it; else
        the estimator default. estimator(1) is forced to poolmass (its h1
        model is h0).
    */
    if "`constant'" != ""         local useconstant = 1
    else if "`exact'" != ""       local useconstant = 0
    else if "`_e_constant'" != "" local useconstant = `_e_constant'
    else                          local useconstant = inlist(`estimator',1,2)

    if "`poolmass'" != ""         local nosplit = 1
    else if "`splitmass'" != ""   local nosplit = 0
    else if "`_e_nosplit'" != ""  local nosplit = `_e_nosplit'
    else                          local nosplit = inlist(`estimator',1,2)
    if `estimator' == 1 local nosplit = 1

    /*
        Estimators 0 and 3 are consistent under exact + splitmass: their
        only biases are the two toggleable ones. With both at their
        default there is nothing to report and asking for it is almost
        certainly a mistake, so error rather than silently return zero.
    */
    if inlist(`estimator',0,3) & !`useconstant' & !`nosplit' {
        di as error "polbunchbias: estimator(`estimator') is consistent under exact integration and a split mass calc -- zero bias by construction. Specify constant (constant-density bias), poolmass (overhang bias), or both; or use estimator(1), estimator(2), or estimator(4)."
        exit 198
    }

    if "`bw'" == "" local bw = 1

    // ----------------------------------------------------------------
    // Common: run the bias calculation
    // ----------------------------------------------------------------

    local ehat = `elasticity'
    local lhat = `lambda'
    local ecur = `elasticity'
    local lcur = `lambda'
    local iter = 0
    local conv = .

    tempname out

    if "`iterate'" != "" {
        local conv = 0
        forvalues k = 1/`maxiter' {
            mata: st_matrix("`out'", pb_fobias_core( ///
                `estimator', `bmodel', `islog_val', `useconstant', `nosplit', ///
                `zstar', `t0', `t1', `lcur', `ecur', ///
                `zlo', `zhi', `zl', `zh', `bw' ///
            ))

            scalar __bias_e = el(`out',1,24)
            scalar __bias_l = el(`out',1,26)

            local enew = `ehat' - scalar(__bias_e)
            local lnew = `lhat' - scalar(__bias_l)
            local diff = max(abs(`enew' - `ecur'), abs(`lnew' - `lcur'))

            /*
                Plain fixed-point substitution (ecur = enew) oscillates
                with growing amplitude for this bias map -- it's a stable
                mapping in the sense that the correct fixed point exists
                and the single-step correction is accurate, but the raw
                iteration's local derivative exceeds 1 in magnitude, so
                naive substitution diverges rather than converging no
                matter how many iterations run. Under-relaxation (moving
                only underrelax() of the way to the candidate each step)
                damps that below the divergence threshold.
            */
            local ecur = `ecur' + `underrelax' * (`enew' - `ecur')
            local lcur = `lcur' + `underrelax' * (`lnew' - `lcur')
            local iter = `k'

            if `diff' < `tolerance' {
                local conv = 1
                continue, break
            }
        }
    }

    mata: st_matrix("`out'", pb_fobias_core( ///
        `estimator', `bmodel', `islog_val', `useconstant', `nosplit', ///
        `zstar', `t0', `t1', `lcur', `ecur', ///
        `zlo', `zhi', `zl', `zh', `bw' ///
    ))

    local names estimator bmodel islog zstar t0 t1 tau lambda elasticity x rho Delta ///
        zlo zhi zL zH dL dR B bias_h bias_B bias_response bias_shift ///
        bias_elasticity bias_slope bias_lambda

    matrix colnames `out' = `names'

    tempname rb
    matrix `rb' = ( ///
        el(`out',1,20), ///
        el(`out',1,25), ///
        el(`out',1,26), ///
        el(`out',1,21), ///
        el(`out',1,22), ///
        el(`out',1,23), ///
        el(`out',1,24) ///
    )
    matrix colnames `rb' = h slope relative_slope number_bunchers marginal_response shift elasticity

    // Display using `rb' directly -- r(b) inside an rclass program still
    // points to the caller's return space, so reading r(b) here would give
    // stale results from whatever command ran before polbunchbias.
    di as text _newline "Polynomial bunching bias estimates: Estimator `estimator'"
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

    // Post return values after display so r() is clean on exit.
    forvalues j = 1/26 {
        local nm : word `j' of `names'
        tempname s`j'
        scalar `s`j'' = el(`out',1,`j')
        return scalar `nm' = `s`j''
    }
    return matrix b = `rb'
    return scalar constant = `useconstant'
    return scalar nosplit = `nosplit'
    return scalar input_elasticity = `ehat'
    return scalar input_lambda = `lhat'
    return scalar corrected_elasticity = `ecur'
    return scalar corrected_lambda = `lcur'
    return scalar iterations = `iter'
    return scalar converged = `conv'

end

capture mata: mata drop pb_fobias_core()
capture mata: mata drop pb_solve_quad()
mata:

real rowvector pb_intP(real scalar lo, real scalar hi)
{
    return((hi-lo, (hi^2-lo^2)/2))
}

real matrix pb_intPP(real scalar lo, real scalar hi)
{
    return((
        hi-lo,              (hi^2-lo^2)/2 \
        (hi^2-lo^2)/2,      (hi^3-lo^3)/3
    ))
}

real scalar pb_solve_quad(real scalar aa, real scalar mm, real scalar S)
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

real rowvector pb_fobias_core(
    real scalar estimator,
    real scalar bmodel,
    real scalar islog,
    real scalar useconstant,
    real scalar nosplit,
    real scalar zstar,
    real scalar t0,
    real scalar t1,
    real scalar lambda,
    real scalar elast,
    real scalar zlo,
    real scalar zhi,
    real scalar zL,
    real scalar zH,
    real scalar bw
)
{
    real scalar tau, Ltau, x, rho, Delta, a, m, r, B
    real scalar dL, dR
    real scalar bias_h, bias_B, bias_resp, bias_shift, bias_e
    real scalar bias_slope, bias_lambda
    real scalar atilde, mtilde, Btilde, rtilde
    real scalar r_const_bias
    real matrix M, GtG, Gtu, biaspar
    real rowvector R
    real scalar lo, hi, L, H, A0, A1, A2
    real scalar gc0, gc1, gn0, gn1
	real scalar q, dqdD
	real scalar gu0, gu1
	real scalar uM, Sbar, Sright
	real scalar lo_right, lo_true, hi_true, trueRightMass, trueMass
	real scalar edge_ovh, overhang
	real scalar delc, dstep, nlsit
	real rowvector beta0
	real rowvector Rall, Rlo, Rhi, Rbar, Jm
	
			real scalar sleft, sright, sright0, hright0
			real scalar m_saez, a_saez


    tau   = (1-t0)/(1-t1)
    Ltau  = ln(tau)
    x     = tau^elast
    rho   = ln(x)
    Delta = x - 1

    if (bw <= 0) bw = 1

    a = 1
    /*
        lambda (the argument) is the EARNINGS elasticity of h0 at the
        cutoff, z* h0'(z*)/h0(z*) -- the scale-free, z*-independent
        parameterisation the caller's relslope() maps to
        (lambda = relslope * zstar, relslope = m/a per unit earnings).

        m below is the slope of h0 in the RUNNING-variable coordinate,
        which differs between the two earnings measures:
          levels : running var is z,     m = (dh0/dz)/a        = lambda/zstar
          logs   : running var is ln z,  m = (dh0/d ln z)/a    = lambda
        Every formula downstream pairs m with running-coordinate
        distances (rho, H, L, hi, lo, ...), so it must be the
        running-coordinate slope. Feeding the per-earnings value into the
        log branch (the old behaviour) made a given relslope() denote a
        density z* times flatter in logs than in levels.
    */
    m = islog ? lambda : lambda / zstar

    lo = zlo - zstar
    hi = zhi - zstar
    L  = zL  - zstar
    H  = zH  - zstar

    if (islog) {
        r = rho
        B = (a*rho + 0.5*m*rho^2)/bw
    }
    else {
        r = zstar*Delta
        B = (a*r + 0.5*m*r^2)/bw
    }

    dL = zstar - zL + (zL-zlo)/2
    dR = zH - zstar + (zhi-zH)/2

    bias_h = bias_B = bias_resp = bias_shift = bias_e = .
    bias_slope = bias_lambda = .

      if (estimator == 4) {
        /*
            Estimator 4: Saez three-region trapezoid approximation.

            Production estimator:
                hminus = mean count in left reference region
                hplus  = mean count in right reference region
                Bsaez  = Hstar - a0*hminus - a1*hplus

            Then saez_transform solves:
                Bsaez = zstar/(2*bw) * (xhat - 1) *
                        (hminus + hplus/xhat)

            This block mirrors that estimator directly instead of using
            the previous dL/dR shortcut.
        */

        real scalar a0s, a1s
        real scalar hminus, hplus, Hstar, Bsaez
        real scalar A, qsaez, disc, xhat, dlogzhat
        real scalar B0

        a0s = -L / bw
        a1s =  H / bw

        if (islog) {
            /*
                Log case:
                    At lambda = 0:
                        hminus = hplus = a
                        B0 = a*rho/bw

                    Right-side log displacement contributes m*rho
                    to the right-side density.
            */

            hminus = a + m*((lo + L)/2)
            hplus  = a + m*((H + hi)/2 + rho)

            Hstar = (a*(H-L) + 0.5*m*(H^2-L^2))/bw + ///
                    (a*rho + 0.5*m*rho^2)/bw + ///
                    m*H*rho/bw

            Bsaez = Hstar - a0s*hminus - a1s*hplus

			/*
				Log case:
					h1(s) = h0(s + rho)

				So hplus corresponds to the counterfactual point:
					(sright + rho, hplus)
			*/


			sleft  = (lo + L)/2
			sright = (H  + hi)/2

			sright0 = sright + rho
			hright0 = hplus

			m_saez = (hright0 - hminus)/(sright0 - sleft)
			a_saez = hminus - m_saez*sleft

			bias_h      = a_saez - a
			bias_slope  = m_saez - m
			bias_lambda = (islog ? 1 : zstar)*(m_saez/a_saez - m/a)

			if (hminus > 0 & hplus > 0) {
				/*
					Constant-density approximation: use hminus alone as
					the counterfactual reference, dropping hplus -- Bsaez
					itself is untouched. At hminus==hplus this coincides
					exactly with the exact formula below (no approximation
					error), same as in saez_transform.
				*/
				if (useconstant) dlogzhat = Bsaez*bw/hminus
				else             dlogzhat = 2*Bsaez*bw/(hminus + hplus)
				xhat     = exp(dlogzhat)
				bias_B    = Bsaez - B
				bias_resp = dlogzhat - rho
				bias_shift = .
				bias_e    = bias_resp / Ltau
			}
        }
        else {
            /*
                Level case.

                Left reference region:
                    hminus = average h0 over [lo,L].

                Right reference region:
                    hplus = average post-reform right density over [H,hi].

                The right-side density under the isoelastic model is,
                to first order in the linear slope,

                    h1(s) = x*a
                            + x*m*(x-1)*zstar
                            + x^2*m*s,

                where s = z - zstar.
            */

            hminus = a + m*((lo + L)/2)

            hplus = x*a + x*m*(x-1)*zstar + ///
                    x^2*m*((H + hi)/2)

			Hstar = (a*(0-L) + 0.5*m*(0^2-L^2))/bw + ///
					(a*r + 0.5*m*r^2)/bw + ///
					(x*H*(a + m*zstar*(x-1)) + 0.5*m*x^2*H^2)/bw
		
            Bsaez = Hstar - a0s*hminus - a1s*hplus

            /*
			Saez implicit counterfactual line.

			hminus is a left-side counterfactual point:
				(sleft, hminus)

			hplus is observed post-bunching right-side density.  To recover
			the corresponding counterfactual h0 point, shift it back:

				observed sright maps to pre-bunching
					sright0 = (x - 1)*zstar + x*sright

				and density rescales by the Jacobian:
					hright0 = hplus/x
		*/

		sleft  = (lo + L)/2
		sright = (H  + hi)/2

		sright0 = (x - 1)*zstar + x*sright
		hright0 = hplus/x

		m_saez = (hright0 - hminus)/(sright0 - sleft)
		a_saez = hminus - m_saez*sleft

		bias_h      = a_saez - a
		bias_slope  = m_saez - m
		bias_lambda = (islog ? 1 : zstar)*(m_saez/a_saez - m/a)

		if (hminus > 0) {
			A = 2*Bsaez*bw/zstar
			/*
				Constant-density approximation: r_const = Bsaez*bw/hminus,
				i.e. xhat = 1 + A/(2*hminus), dropping hplus (and the
				quadratic in xhat) entirely -- coincides with the exact
				solve below only to first order in the response, even at
				hminus==hplus, mirroring saez_transform.
			*/
			if (useconstant) {
				xhat = 1 + A/(2*hminus)
				bias_B    = Bsaez - B
				bias_resp = zstar*(xhat - x)
				bias_shift = xhat - x
				bias_e    = ln(xhat)/Ltau - elast
			}
			else {
				qsaez = hplus - hminus - A
				disc  = qsaez^2 + 4*hminus*hplus
				if (disc >= 0) {
					xhat = (-qsaez + sqrt(disc))/(2*hminus)
					if (xhat > 0) {
						bias_B    = Bsaez - B
						bias_resp = zstar*(xhat - x)
						bias_shift = xhat - x
						bias_e    = ln(xhat)/Ltau - elast
					}
				}
			}
		}
		}
		
        return((estimator,bmodel,islog,zstar,t0,t1,tau,lambda,elast,x,rho,Delta, ///
            zlo,zhi,zL,zH,dL,dR,B,bias_h,bias_B,bias_resp,bias_shift,bias_e, ///
            bias_slope,bias_lambda))
    }
    if (estimator == 1) {
        M = pb_intPP(lo,L) + pb_intPP(H,hi)

        if (islog) {
            Gtu = pb_intP(H,hi)' * (m*rho)
        }
        else {
            A0 = hi-H
            A1 = (hi^2-H^2)/2
            A2 = (hi^3-H^3)/3

            gn0 = a*(x-1) + x*m*(x-1)*zstar
            gn1 = m*(x^2-1)

            Gtu = J(2,1,0)
            Gtu[1] = gn0*A0 + gn1*A1
            Gtu[2] = gn0*A1 + gn1*A2
        }

        biaspar = invsym(M)*Gtu

        bias_h = biaspar[1]
        bias_slope = biaspar[2]
        // prefactor turns d(m/a) into d(lambda): lambda = zstar*(m/a) in
        // levels, but lambda = (m/a) directly in logs (running var already ln z)
        bias_lambda = (islog ? 1 : zstar) * (bias_slope/a - (m/a)*(bias_h/a))

        R = pb_intP(L,H)/bw
        bias_B = -(R*biaspar)[1]

        /*
            Overhang / compression term. The naive Bhat is (observed mass
            in E) - (fitted h0 integrated over E). The observed mass in E
            is NOT the counterfactual mass over E: the behavioural response
            compresses counterfactual mass from ABOVE zH into the observed
            window, out to the stretched upper edge s0 = x*H + r (levels)
            or H + rho (logs) -- the same edge Chetty's mass row uses
            (hi_true / Rall below). So even with a perfect counterfactual
            fit, Bhat*bw = int_{H}^{edge} h0, not the true int_0^{r} h0 =
            B*bw. -(R*biaspar) above is only the fit-error piece; add the
            remaining [ int_{H}^{edge} h0 ]/bw - B. The size depends on the
            UPPER edge only (zH vs z*), not on zL: it is present whenever
            the excluded window reaches above the cutoff (i.e. always in
            practice). It nearly cancels for a narrow, steeply sloped level
            window (why the level baseline looked exact) but is first-order
            in rho for logs and grows with the fitting window in levels.
            Guarded on L<0<H -- the naive Bhat is only meaningful when the
            excluded window straddles the cutoff.

            Not gated on nosplit: the naive estimator has no h1 model
            distinct from h0, so "split the window at the kink" and "pool
            h0 over the whole window" are the same operation -- poolmass is
            definitionally on for estimator(1), splitmass a no-op (the
            .ado forces nosplit=1 for estimator 1).
        */
        if (L < 0 & H > 0) {
            if (islog) edge_ovh = H + rho
            else       edge_ovh = x*H + r
            overhang = (a*(edge_ovh - H) + 0.5*m*(edge_ovh^2 - H^2))/bw - B
            bias_B = bias_B + overhang
        }

        atilde = a + bias_h
        mtilde = m + bias_slope
        Btilde = B + bias_B

        if (useconstant) {
            /*
                Constant-density mode: the estimator recovers the response
                as rhat = Bhat*bw / h0hat(z*) = Btilde*bw/atilde, dropping
                the fitted slope over the bunching region (polbunch.ado's
                MR = B*bw/m with m = h0 at the cutoff). Use that exact ratio
                and the same r->e conversion as the exact branch below -- an
                earlier version linearized this as Delta*(bias_B/B -
                bias_h/a), which is a further approximation on top of the
                constant-density one and left the constant cells flagged.
            */
            rtilde = Btilde*bw/atilde
        }
        else {
            rtilde = pb_solve_quad(atilde, mtilde, Btilde*bw)
        }

        bias_resp = rtilde - r
        if (islog) {
            bias_shift = .
            bias_e = bias_resp/Ltau
        }
        else {
            bias_shift = rtilde/zstar - Delta
            bias_e = ln(1 + rtilde/zstar)/Ltau - elast
        }

        return((estimator,bmodel,islog,zstar,t0,t1,tau,lambda,elast,x,rho,Delta, ///
            zlo,zhi,zL,zH,dL,dR,B,bias_h,bias_B,bias_resp,bias_shift,bias_e, ///
            bias_slope,bias_lambda))
    }

 if (estimator == 2) {
    if (islog) {
        /*
            Estimator 2, log-earnings case.

            Chetty's own production restriction is unit-agnostic -- it
            says h1=h0/(1+Delta) at the same argument s, same q=1/(1+Delta)
            scaling, regardless of whether s is a level or log deviation --
            so the shape-fitting design (GtG's q, dqdD structure below) is
            identical to the level case. What changes is the truth being
            fit against: under the log-earnings iso-elastic model, the
            true right-side density is the pure translation h0(s+rho), not
            the proportional relocation x*h0(xs) of the level case, and the
            true excluded-region mass is the integral of h0 extended by
            rho past H -- bunchers collapse exactly onto the cutoff,
            refilling the window from what would otherwise have sat just
            above it -- rather than the level mass row's separate
            Delta*int_0^zbar h0 term standing in for the excess mass.
        */
        beta0 = (a, m)

        A0 = hi-H
        A1 = (hi^2-H^2)/2
        A2 = (hi^3-H^3)/3

        /*
            gn0/gn1: true right density minus PLAIN h0. In logs the
            relocation is the pure translation h0(s+rho), so true h1 - h0 =
            m*rho (level shift only, no slope jump). The estimator-2
            residual "true h1 - q*h0" adds (1-q)h0, with q = 1/(1+delta) at
            the NLS pseudo-true delta from the iteration below.
        */
        gn0 = m*rho
        gn1 = 0

        Rall = pb_intP(L, H+rho)/bw

        Rlo = J(1,2,0)
        if (L < 0) {
            if (H <= 0) Rlo = pb_intP(L,H)/bw
            else        Rlo = pb_intP(L,0)/bw
        }
        Rhi = J(1,2,0)
        if (H > 0) {
            if (L >= 0) Rhi = pb_intP(L,H)/bw
            else        Rhi = pb_intP(0,H)/bw
        }
        Rbar   = pb_intP(0,hi)/bw
        Sbar   = Rbar * beta0'
        Sright = Rhi  * beta0'

        /*
            True content of [L,H] in logs: bunchers collapse onto the
            cutoff, so the window holds the counterfactual mass out to
            H+rho -- no separate B addend (unlike the level case).
        */
        trueMass = Rall*beta0'

        /* NLS fixed point for delta (see the level branch). */
        delc = Delta
        for (nlsit = 1; nlsit <= 80; nlsit++) {
            q    = 1/(1+delc)
            dqdD = -1/((1+delc)^2)

            gu0 = gn0 + (1-q)*a
            gu1 = gn1 + (1-q)*m

            GtG = J(3,3,0)
            Gtu = J(3,1,0)

            GtG[1..2,1..2] = pb_intPP(lo,L)

            GtG[1,1] = GtG[1,1] + q^2*A0
            GtG[1,2] = GtG[1,2] + q^2*A1
            GtG[2,1] = GtG[1,2]
            GtG[2,2] = GtG[2,2] + q^2*A2

            GtG[1,3] = q*dqdD*(a*A0 + m*A1)
            GtG[3,1] = GtG[1,3]
            GtG[2,3] = q*dqdD*(a*A1 + m*A2)
            GtG[3,2] = GtG[2,3]
            GtG[3,3] = dqdD^2 * (a^2*A0 + 2*a*m*A1 + m^2*A2)

            Gtu[1] = q*(gu0*A0 + gu1*A1)
            Gtu[2] = q*(gu0*A1 + gu1*A2)
            Gtu[3] = dqdD*(a*gu0*A0 + (a*gu1 + m*gu0)*A1 + m*gu1*A2)

            Jm = J(1,3,0)
            Jm[1,1..2] = Rlo + q*Rhi + delc*Rbar
            Jm[1,3]    = dqdD*Sright + Sbar

            uM = trueMass - ((Rlo + q*Rhi + delc*Rbar)*beta0')

            GtG = GtG + Jm'Jm
            Gtu = Gtu + Jm'*uM

            biaspar = invsym(GtG)*Gtu

            dstep = 0.7 * biaspar[3]
            delc  = delc + dstep
            if (1 + delc <= 1e-6) delc = -1 + 1e-6
            if (abs(dstep) < 1e-11) break
        }

        bias_h = biaspar[1]
        bias_slope = biaspar[2]
        // prefactor turns d(m/a) into d(lambda): lambda = zstar*(m/a) in
        // levels, but lambda = (m/a) directly in logs (running var already ln z)
        bias_lambda = (islog ? 1 : zstar) * (bias_slope/a - (m/a)*(bias_h/a))

        /*
            Reduced-form Bhat = Hstar - int_E h0hat: step-1 fit error plus
            the step-2 term of Section~sec:overhang. poolmass subtracts
            h0hat both sides (overhang, same as naive); splitmass subtracts
            h0hat/x on the right, replacing the overhang with
            int_{z*}^{zH}(h1_true - h0/x) = overhang + (1-1/x) int_0^H h0.
        */
        R = pb_intP(L,H)/bw
        bias_B = -(R*biaspar[1..2])[1]

        if (L < 0 & H > 0) {
            edge_ovh = H + rho
            overhang = (a*(edge_ovh - H) + 0.5*m*(edge_ovh^2 - H^2))/bw - B
            if (nosplit) bias_B = bias_B + overhang
            else           bias_B = bias_B + overhang + (1 - 1/x)*(a*H + 0.5*m*H^2)/bw
        }

        atilde = a + bias_h
        mtilde = m + bias_slope
        Btilde = B + bias_B

        if (useconstant) rtilde = Btilde*bw/atilde
        else             rtilde = pb_solve_quad(atilde, mtilde, Btilde*bw)

        bias_resp  = rtilde - rho
        bias_shift = .
        bias_e     = bias_resp/Ltau

        return((estimator,bmodel,islog,zstar,t0,t1,tau,lambda,elast,x,rho,Delta, ///
            zlo,zhi,zL,zH,dL,dR,B,bias_h,bias_B,bias_resp,bias_shift,bias_e, ///
            bias_slope,bias_lambda))
    }

		/*
			Estimator 2, level case.

			Production restriction:
				h1 = h0 / (1 + delta)

			But production estimation also includes the Chetty mass row:
				Hstar = counterfactual mass in excluded region
						+ delta * int_{zstar}^{zbar} h0(z) dz / bw

			With bmodel(0), the final reported B is still reduced-form:
				B_RF = Hstar - int_{zL}^{zH} h0(z) dz / bw.
		*/


		beta0 = (a, m)

		A0 = hi-H
		A1 = (hi^2-H^2)/2
		A2 = (hi^3-H^3)/3

		/*
			True right-side density minus PLAIN h0 (the estimator-1
			residual). The estimator-2 residual "true h1 - q*h0" adds
			(1-q)*h0, with q = 1/(1+delta) at the NLS pseudo-true delta
			found by the iteration below -- NOT q = 1/x. polbunch's
			estimator 2 profiles delta jointly; it converges to a
			window-and-zbar-dependent value (near 0 for wide support,
			~1.5-2*Delta for narrow), not to Delta.
		*/
		gn0 = a*(x-1) + x*m*(x-1)*zstar
		gn1 = m*(x^2-1)

		/*
			delta-independent window integrals for the Chetty mass row.
			Rbar = int_{z*}^{zbar} P_K / bw, zbar = zhi (the compressed
			observed fitting bound, as passed by the caller).
		*/
		Rlo = J(1,2,0)
		if (L < 0) {
			if (H <= 0) Rlo = pb_intP(L,H)/bw
			else        Rlo = pb_intP(L,0)/bw
		}
		Rhi = J(1,2,0)
		if (H > 0) {
			if (L >= 0) Rhi = pb_intP(L,H)/bw
			else        Rhi = pb_intP(0,H)/bw
		}
		Rbar   = pb_intP(0,hi)/bw
		Sbar   = Rbar * beta0'
		Sright = Rhi  * beta0'

		/*
			True content of [L,H] under the exact iso-elastic relocation
			s0 = x*s + r (DGP; delta-independent). Below 0 untouched
			(through Rlo); point mass B at the cutoff; above 0 the mass
			over (max(L,0),H] maps from counterfactual (x*max(L,0)+r, x*H+r].
		*/
		lo_right = max((L,0))
		if (H > lo_right) {
			lo_true = x*lo_right + r
			hi_true = x*H + r
			trueRightMass = (a*(hi_true-lo_true) + 0.5*m*(hi_true^2-lo_true^2))/bw
		}
		else trueRightMass = 0
		trueMass = Rlo*beta0' + trueRightMass
		if (L < 0 & H > 0) trueMass = trueMass + B

		/*
			NLS fixed point (equivalent to Chetty's iterative integration
			constraint). Linearise the beta-block around beta_true and the
			delta-block around the current delc; the Newton step's delta
			component drives delc -> delta*. At convergence biaspar[3] ~ 0
			and biaspar[1..2] = (b_h, b_slope). Everything downstream (Bhat,
			response, elasticity) uses only b_h, b_slope -- delta* is a
			nuisance that never propagates.
		*/
		delc = Delta
		for (nlsit = 1; nlsit <= 80; nlsit++) {
			q    = 1/(1+delc)
			dqdD = -1/((1+delc)^2)

			gu0 = gn0 + (1-q)*a
			gu1 = gn1 + (1-q)*m

			GtG = J(3,3,0)
			Gtu = J(3,1,0)

			GtG[1..2,1..2] = pb_intPP(lo,L)

			GtG[1,1] = GtG[1,1] + q^2*A0
			GtG[1,2] = GtG[1,2] + q^2*A1
			GtG[2,1] = GtG[1,2]
			GtG[2,2] = GtG[2,2] + q^2*A2

			GtG[1,3] = q*dqdD*(a*A0 + m*A1)
			GtG[3,1] = GtG[1,3]
			GtG[2,3] = q*dqdD*(a*A1 + m*A2)
			GtG[3,2] = GtG[2,3]
			GtG[3,3] = dqdD^2 * (a^2*A0 + 2*a*m*A1 + m^2*A2)

			Gtu[1] = q*(gu0*A0 + gu1*A1)
			Gtu[2] = q*(gu0*A1 + gu1*A2)
			Gtu[3] = dqdD*(a*gu0*A0 + (a*gu1 + m*gu0)*A1 + m*gu1*A2)

			Jm = J(1,3,0)
			Jm[1,1..2] = Rlo + q*Rhi + delc*Rbar
			Jm[1,3]    = dqdD*Sright + Sbar

			uM = trueMass - ((Rlo + q*Rhi + delc*Rbar)*beta0')

			GtG = GtG + Jm'Jm
			Gtu = Gtu + Jm'*uM

			biaspar = invsym(GtG)*Gtu

			dstep = 0.7 * biaspar[3]
			delc  = delc + dstep
			if (1 + delc <= 1e-6) delc = -1 + 1e-6
			if (abs(dstep) < 1e-11) break
		}

		bias_h = biaspar[1]
		bias_slope = biaspar[2]
		bias_lambda = (islog ? 1 : zstar) * (bias_slope/a - (m/a)*(bias_h/a))

		if (bmodel == 1) {
			/*
				Final reported B/response is model-implied: the response is
				the NLS pseudo-true delta itself, so the shift bias is
				delta* - Delta (biaspar[3] is ~0 at the fixed point).
			*/
			bias_shift = delc - Delta
			bias_resp  = zstar*bias_shift
			bias_B     = .
			bias_e     = ln(1 + Delta + bias_shift)/Ltau - elast
		}
		else {
			/*
				Reduced-form Bhat = Hstar - int_E h0hat: step-1 fit error
				plus the step-2 term of Section~sec:overhang. poolmass
				subtracts h0hat on both sides (the overhang, same term as
				the naive estimator); splitmass subtracts h0hat/x on the
				right (Chetty's own h1 model), replacing the overhang with
				int_{z*}^{zH}(h1_true - h0/x) = overhang + (1-1/x) int_0^H h0.
			*/
			R = pb_intP(L,H)/bw
			bias_B = -(R*biaspar[1..2])[1]

			if (L < 0 & H > 0) {
				edge_ovh = x*H + r
				overhang = (a*(edge_ovh - H) + 0.5*m*(edge_ovh^2 - H^2))/bw - B
				if (nosplit) bias_B = bias_B + overhang
				else           bias_B = bias_B + overhang + (1 - 1/x)*(a*H + 0.5*m*H^2)/bw
			}

			atilde = a + bias_h
			mtilde = m + bias_slope
			Btilde = B + bias_B

			if (useconstant) rtilde = Btilde*bw/atilde
			else             rtilde = pb_solve_quad(atilde, mtilde, Btilde*bw)

			bias_resp  = rtilde - r
			bias_shift = rtilde/zstar - Delta
			bias_e     = ln(1 + rtilde/zstar)/Ltau - elast
		}

		return((estimator,bmodel,islog,zstar,t0,t1,tau,lambda,elast,x,rho,Delta, ///
			zlo,zhi,zL,zH,dL,dR,B,bias_h,bias_B,bias_resp,bias_shift,bias_e, ///
			bias_slope,bias_lambda))
	}

    if (estimator == 3 | estimator == 0) {
        /*
            Estimators 3 and 0 have no fitting-stage bias: (a,m) are the
            truth -- estimator 3 because it imposes the correct
            shift-and-stretch restriction, estimator 0 because it fits h0
            from the left data only (uncontaminated) and h1 freely from the
            right data. The only biases they carry are the two downstream
            researcher choices, toggled independently of the estimator:

              nosplit  (axis 2) -- read Bhat off as raw excess over h0
                          pooled across the whole window instead of
                          splitting at the kink. Adds the same overhang
                          term as the naive estimator, at the (here exact)
                          h0.
              useconstant (axis 3) -- recover the response as Bhat*bw/h0(z*)
                          instead of inverting the integral.

            With both off the branch returns zero (the .ado guards against
            reaching here in that case). With both on it returns
            beta_ovh + beta_const, plus the second-order cross term from
            compounding both through the nonlinear level r->e map (in logs
            that map is linear, so there the two are exactly additive).
        */
        bias_h = bias_slope = bias_lambda = 0

        bias_B = 0
        if (nosplit & L < 0 & H > 0) {
            if (islog) edge_ovh = H + rho
            else       edge_ovh = x*H + r
            bias_B = (a*(edge_ovh - H) + 0.5*m*(edge_ovh^2 - H^2))/bw - B
        }

        Btilde = B + bias_B
        if (useconstant) {
            /*
                rhat_const = Bhat*bw / h0hat(z*) = Btilde*bw/a (h0hat(z*) = a
                exactly for estimator 3). In the RESPONSE, axes 2 and 3 are
                exactly additive: rhat_const - r = 0.5*(m/a)*r^2 (constant-
                approx) + Btilde_ovh*bw/a (overhang). The level r->e step is
                nonlinear, so beta_e picks up a second-order cross term;
                logs stay additive.
            */
            rtilde = Btilde*bw/a
        }
        else {
            rtilde = pb_solve_quad(a, m, Btilde*bw)
        }

        if (bias_B == 0 & !useconstant) {
            bias_resp = bias_shift = bias_e = 0
        }
        else {
            bias_resp = rtilde - r
            if (islog) {
                bias_shift = .
                bias_e     = bias_resp/Ltau
            }
            else {
                bias_shift = rtilde/zstar - Delta
                bias_e     = ln(1 + rtilde/zstar)/Ltau - elast
            }
        }

        return((estimator,bmodel,islog,zstar,t0,t1,tau,lambda,elast,x,rho,Delta, ///
            zlo,zhi,zL,zH,dL,dR,B,bias_h,bias_B,bias_resp,bias_shift,bias_e, ///
            bias_slope,bias_lambda))
    }

    return(J(1,26,.))
}

end