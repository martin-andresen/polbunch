capture program drop polbunchbias
program define polbunchbias, rclass
    version 16.0

    syntax [, ESTimator(numlist max=1) ZSTAR(numlist max=1) ///
        T0(numlist max=1) T1(numlist max=1) ///
        LAMBDA(numlist max=1) ELasticity(numlist max=1) ///
        ZLO(numlist max=1) ZHI(numlist max=1) ///
        ZL(numlist max=1) ZH(numlist max=1) ///
        BMODEL(integer 0) LOG ///
        BW(numlist max=1) ITERate TOLerance(real 1e-10) ///
        MAXITER(integer 100) CONstant ]

    // ----------------------------------------------------------------
    // Mode detection: count how many of the 9 primary required
    // options were supplied.  lambda is handled separately.
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

        // lambda: compute from e(b) unless overridden by user
        if "`lambda'" == "" {
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
                        local _hval  = `_hval'  + `_bk' * `_cest'^`_k'
                        local _dhval = `_dhval' + `_k' * `_bk' * `_cest'^(`_k'-1)
                    }
                }
                if `_hval' > 0 {
                    local lambda = `_dhval' * `zstar' / (`_xscale' * `_hval')
                }
                else local lambda = 0
            }
        }

        if "`bw'" == "" local bw = `_bw'
    }
    else if `ncore' == 9 {
        // ============================================================
        // EXPLICIT MODE: all 9 primary options supplied
        // ============================================================
        if "`lambda'"     == "" local lambda = 0
        if "`bw'"         == "" local bw = 1
        local islog_val = ("`log'" != "")
    }
    else {
        // Partial specification — error
        local missing_opts
        foreach v in estimator zstar t0 t1 elasticity zlo zhi zl zh {
            if "``v''" == "" local missing_opts `missing_opts' `v'()
        }
        di as error "Specify all required options, or none to use polbunch results from e()"
        di as error "Missing: `missing_opts'"
        exit 198
    }

    if "`bw'" == "" local bw = 1

    // ----------------------------------------------------------------
    // Common: run the bias calculation
    // ----------------------------------------------------------------
    local useconstant = ("`constant'" != "")

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
                `estimator', `bmodel', `islog_val', `useconstant', ///
                `zstar', `t0', `t1', `lcur', `ecur', ///
                `zlo', `zhi', `zl', `zh', `bw' ///
            ))

            scalar __bias_e = el(`out',1,24)
            scalar __bias_l = el(`out',1,26)

            local enew = `ehat' - scalar(__bias_e)
            local lnew = `lhat' - scalar(__bias_l)
            local diff = max(abs(`enew' - `ecur'), abs(`lnew' - `lcur'))

            local ecur = `enew'
            local lcur = `lnew'
            local iter = `k'

            if `diff' < `tolerance' {
                local conv = 1
                continue, break
            }
        }
    }

    mata: st_matrix("`out'", pb_fobias_core( ///
        `estimator', `bmodel', `islog_val', `useconstant', ///
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
    real matrix M, GtG, Gtu, biaspar
    real rowvector R
    real scalar lo, hi, L, H, A0, A1, A2
    real scalar gc0, gc1, gn0, gn1
	real scalar q, dqdD
	real scalar gu0, gu1
	real scalar uM, Sbar, Sright
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
    m = lambda / zstar

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
			bias_lambda = zstar*(m_saez/a_saez - m/a)

			if (hminus > 0 & hplus > 0) {
				dlogzhat = 2*Bsaez*bw/(hminus + hplus)
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
		bias_lambda = zstar*(m_saez/a_saez - m/a)

		if (hminus > 0) {
			A     = 2*Bsaez*bw/zstar
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
        bias_lambda = zstar * (bias_slope/a - (m/a)*(bias_h/a))

        R = pb_intP(L,H)/bw
        bias_B = -(R*biaspar)[1]

        if (useconstant) {
            if (islog) {
                bias_resp = rho*(bias_B/B - bias_h/a)
                bias_shift = .
                bias_e = bias_resp/Ltau
            }
            else {
                bias_shift = Delta*(bias_B/B - bias_h/a)
                bias_resp  = zstar*bias_shift
                bias_e     = bias_shift/(x*Ltau)
            }
        }
        else {
            atilde = a + bias_h
            mtilde = m + bias_slope
            Btilde = B + bias_B
            rtilde = pb_solve_quad(atilde, mtilde, Btilde*bw)

            bias_resp = rtilde - r
            if (islog) {
                bias_shift = .
                bias_e = bias_resp/Ltau
            }
            else {
                bias_shift = rtilde/zstar - Delta
                bias_e = ln(1 + rtilde/zstar)/Ltau - elast
            }
        }

        return((estimator,bmodel,islog,zstar,t0,t1,tau,lambda,elast,x,rho,Delta, ///
            zlo,zhi,zL,zH,dL,dR,B,bias_h,bias_B,bias_resp,bias_shift,bias_e, ///
            bias_slope,bias_lambda))
    }

 if (estimator == 2) {
    if (islog) {
        M = pb_intPP(lo,L) + pb_intPP(H,hi)
        Gtu = pb_intP(H,hi)' * (m*rho)
        biaspar = invsym(M)*Gtu

        bias_h = biaspar[1]
        bias_slope = biaspar[2]
        bias_lambda = zstar * (bias_slope/a - (m/a)*(bias_h/a))

        R = pb_intP(L,H)/bw
        bias_B = -(R*biaspar)[1]

        if (useconstant) {
            bias_resp  = rho*(bias_B/B - bias_h/a)
            bias_shift = .
            bias_e     = bias_resp/Ltau
        }
        else {
            atilde = a + bias_h
            mtilde = m + bias_slope
            Btilde = B + bias_B
            rtilde = pb_solve_quad(atilde, mtilde, Btilde*bw)

            bias_resp  = rtilde - rho
            bias_shift = .
            bias_e     = bias_resp/Ltau
        }

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


		q     = 1/x
		dqdD  = -1/(x^2)
		beta0 = (a, m)

		A0 = hi-H
		A1 = (hi^2-H^2)/2
		A2 = (hi^3-H^3)/3

		/*
			True right-side density minus estimator-2 model right-side density.

			Existing estimator-1 residual:
				true h1 - h0 = gn0 + gn1*s.

			Estimator-2 residual:
				true h1 - q*h0
				= (true h1 - h0) + (1-q)h0.
		*/

		gn0 = a*(x-1) + x*m*(x-1)*zstar
		gn1 = m*(x^2-1)

		gu0 = gn0 + (1-q)*a
		gu1 = gn1 + (1-q)*m

		GtG = J(3,3,0)
		Gtu = J(3,1,0)

		/*
			Left-side density rows.
		*/
		GtG[1..2,1..2] = GtG[1..2,1..2] + pb_intPP(lo,L)

		/*
			Right-side density rows.
			Design columns are:
				beta columns: q*[1,s]
				delta column: dqdD*h0(s)
		*/
		GtG[1,1] = GtG[1,1] + q^2*A0
		GtG[1,2] = GtG[1,2] + q^2*A1
		GtG[2,1] = GtG[1,2]
		GtG[2,2] = GtG[2,2] + q^2*A2

		GtG[1,3] = GtG[1,3] + q*dqdD*(a*A0 + m*A1)
		GtG[3,1] = GtG[1,3]

		GtG[2,3] = GtG[2,3] + q*dqdD*(a*A1 + m*A2)
		GtG[3,2] = GtG[2,3]

		GtG[3,3] = GtG[3,3] + dqdD^2 * ///
			(a^2*A0 + 2*a*m*A1 + m^2*A2)

		Gtu[1] = q*(gu0*A0 + gu1*A1)
		Gtu[2] = q*(gu0*A1 + gu1*A2)
		Gtu[3] = dqdD * ///
			(a*gu0*A0 + (a*gu1 + m*gu0)*A1 + m*gu1*A2)

		/*
			Add the production estimator-2 mass row.

			Hstar true is represented as:
				int_{zL}^{zH} h0 dz / bw + B_true.

			Model mass row is:
				counterfactual excluded mass under estimator 2
				+ delta * int_{zstar}^{zbar} h0 dz / bw.

			Here zbar is taken to be zhi, matching the usual upper support
			used in this fobias setup.
		*/

		Rall = pb_intP(L,H)/bw

		Rlo = J(1,2,0)
		if (L < 0) {
			if (H <= 0) {
				Rlo = pb_intP(L,H)/bw
			}
			else {
				Rlo = pb_intP(L,0)/bw
			}
		}

		Rhi = J(1,2,0)
		if (H > 0) {
			if (L >= 0) {
				Rhi = pb_intP(L,H)/bw
			}
			else {
				Rhi = pb_intP(0,H)/bw
			}
		}

		Rbar = pb_intP(0,hi)/bw

		Sbar   = Rbar * beta0'
		Sright = Rhi  * beta0'

		/*
			Derivative of the mass row with respect to beta and delta.
		*/
		Jm = J(1,3,0)
		Jm[1,1..2] = Rlo + q*Rhi + Delta*Rbar
		Jm[1,3]    = dqdD*Sright + Sbar

		/*
			Residual in the mass row at the true beta and true Delta.
		*/
		uM = (Rall*beta0' + B) - ((Rlo + q*Rhi + Delta*Rbar)*beta0')

		GtG = GtG + Jm'Jm
		Gtu = Gtu + Jm'*uM

		biaspar = invsym(GtG)*Gtu

		bias_h = biaspar[1]
		bias_slope = biaspar[2]
		bias_lambda = zstar * (bias_slope/a - (m/a)*(bias_h/a))

		if (bmodel == 1) {
			/*
				Final reported B/response is model-implied.
			*/
			bias_shift = biaspar[3]
			bias_resp  = zstar*bias_shift
			bias_B     = .
			bias_e     = ln(1 + Delta + bias_shift)/Ltau - elast
		}
		else {
			/*
				Final reported B is reduced-form:
					B_RF = Hstar - int h0hat over excluded region.
			*/
			R = pb_intP(L,H)/bw
			bias_B = -(R*biaspar[1..2])[1]

			if (useconstant) {
				bias_shift = Delta*(bias_B/B - bias_h/a)
				bias_resp  = zstar*bias_shift
				bias_e     = bias_shift/(x*Ltau)
			}
			else {
				atilde = a + bias_h
				mtilde = m + bias_slope
				Btilde = B + bias_B
				rtilde = pb_solve_quad(atilde, mtilde, Btilde*bw)

				bias_resp  = rtilde - r
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