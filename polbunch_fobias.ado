capture program drop polbunch_fobias
program define polbunch_fobias, rclass
    version 16.0

    syntax , ESTimator(integer) ZSTAR(numlist max=1) ///
        T0(numlist max=1) T1(numlist max=1) ///
        LAMBDA(numlist max=1) ELasticity(numlist max=1) ///
        ZLO(numlist max=1) ZHI(numlist max=1) ///
        ZL(numlist max=1) ZH(numlist max=1) ///
        [ BMODEL(integer 0) LOG DL(numlist max=1) DR(numlist max=1) ///
          BW(numlist max=1) ITERate TOLerance(real 1e-10) ///
          MAXITER(integer 100) CONstant ]

    if "`bw'" == "" local bw = 1
    if "`dl'" == "" local dl = .
    if "`dr'" == "" local dr = .

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
                `estimator', `bmodel', ("`log'"!=""), `useconstant', ///
                `zstar', `t0', `t1', `lcur', `ecur', ///
                `zlo', `zhi', `zl', `zh', `dl', `dr', `bw' ///
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
        `estimator', `bmodel', ("`log'"!=""), `useconstant', ///
        `zstar', `t0', `t1', `lcur', `ecur', ///
        `zlo', `zhi', `zl', `zh', `dl', `dr', `bw' ///
    ))

    local names estimator bmodel islog zstar t0 t1 tau lambda elasticity x rho Delta ///
        zlo zhi zL zH dL dR B bias_h bias_B bias_response bias_shift ///
        bias_elasticity bias_slope bias_lambda

    matrix colnames `out' = `names'

    forvalues j = 1/26 {
        local nm : word `j' of `names'
        tempname s`j'
        scalar `s`j'' = el(`out',1,`j')
        return scalar `nm' = `s`j''
    }

    tempname rb
	matrix `rb' = ( ///
		el(`out',1,19), ///
		el(`out',1,25), ///
		el(`out',1,26), ///
		el(`out',1,20), ///
		el(`out',1,21), ///
		el(`out',1,22), ///
		el(`out',1,24) ///
	)
	matrix colnames `rb' = bias_h bias_slope bias_lambda bias_B ///
		bias_response bias_shift bias_elasticity
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
    real scalar dL,
    real scalar dR,
    real scalar bw
)
{
    real scalar tau, Ltau, x, rho, Delta, a, m, r, B
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

    if (dL >= .) dL = zstar - zL + (zL-zlo)/2
    if (dR >= .) dR = zH - zstar + (zhi-zH)/2

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
                    (a*rho + 0.5*m*rho^2)/bw

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

            Hstar = (a*(H-L) + 0.5*m*(H^2-L^2))/bw + ///
                    (a*r + 0.5*m*r^2)/bw

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