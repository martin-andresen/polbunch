capture program drop bunch_fobias_core
program define bunch_fobias_core, rclass
    version 16.0

	syntax , ESTimator(integer) ZSTAR(numlist max=1) ///
		T0(numlist max=1) T1(numlist max=1) ///
		LAMBDA(numlist max=1) ELasticity(numlist max=1) ///
		ZLO(numlist max=1) ZHI(numlist max=1) ///
		ZL(numlist max=1) ZH(numlist max=1) ///
		[ BMODEL(integer 0) LOG DL(numlist max=1) DR(numlist max=1) ///
		  BW(numlist max=1) ITERate TOLerance(real 1e-10) MAXITER(integer 100) ]
		  
    if "`bw'" == "" local bw = 1
    if "`dl'" == "" local dl = .
    if "`dr'" == "" local dr = .

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
                `estimator', `bmodel', ("`log'"!=""), ///
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
        `estimator', `bmodel', ("`log'"!=""), ///
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
        r(bias_h), ///
        r(bias_slope), ///
        r(bias_lambda), ///
        r(bias_B), ///
        r(bias_response), ///
        r(bias_shift), ///
        r(bias_elasticity) ///
    )
    matrix colnames `rb' = bias_h bias_slope bias_lambda bias_B ///
        bias_response bias_shift bias_elasticity
    return matrix b = `rb'

    return scalar input_elasticity = `ehat'
    return scalar input_lambda = `lhat'
    return scalar corrected_elasticity = `ecur'
    return scalar corrected_lambda = `lcur'
    return scalar iterations = `iter'
    return scalar converged = `conv'

    di as txt "First-order bunching bias"
    di as txt "  estimator             = " as res r(estimator)
    di as txt "  bmodel                = " as res r(bmodel)
    di as txt "  log                   = " as res r(islog)
    di as txt "  x                     = " as res %10.5f r(x)
    di as txt "  relative slope lambda = " as res %10.5f r(lambda)
    di as txt "  bias slope            = " as res %10.5f r(bias_slope)
    di as txt "  bias lambda           = " as res %10.5f r(bias_lambda)
    di as txt "  bias response         = " as res %10.5f r(bias_response)
    di as txt "  bias shift            = " as res %10.5f r(bias_shift)
    di as txt "  bias elasticity       = " as res %10.5f r(bias_elasticity)

    if "`iterate'" != "" {
        di as txt "  input elasticity      = " as res %10.5f r(input_elasticity)
        di as txt "  corrected elasticity  = " as res %10.5f r(corrected_elasticity)
        di as txt "  input lambda          = " as res %10.5f r(input_lambda)
        di as txt "  corrected lambda      = " as res %10.5f r(corrected_lambda)
        di as txt "  iterations            = " as res r(iterations)
        di as txt "  converged             = " as res r(converged)
    }
end

capture mata: mata drop pb_fobias_core()
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

real rowvector pb_fobias_core(
    real scalar estimator,
    real scalar bmodel,
    real scalar islog,
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
    real matrix M, GtG, Gtu, biaspar
    real rowvector R
    real scalar lo, hi, L, H, A0, A1, A2, kg
	real scalar gc0, gc1

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
        if (islog) {
            bias_resp  = -(m/a)*rho*(dR-dL)
            bias_shift = .
            bias_e     = bias_resp/Ltau
        }
        else {
            bias_resp  = -(m/a)*r*(x*dR-dL)
            bias_shift = bias_resp/zstar
            bias_e     = bias_shift/(x*Ltau)
        }

        bias_slope  = 0
        bias_lambda = 0

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

            Gtu = J(2,1,0)
			Gtu[1] = (a*(x-1) + x*m*(x-1)*zstar)*A0 ///
				   + m*(x^2-1)*A1

			Gtu[2] = (a*(x-1) + x*m*(x-1)*zstar)*A1 ///
				   + m*(x^2-1)*A2
        }

        biaspar = invsym(M)*Gtu

        bias_h = biaspar[1]
        bias_slope = biaspar[2]
        bias_lambda = zstar * (bias_slope/a - (m/a)*(bias_h/a))

        R = pb_intP(L,H)/bw
        bias_B = -(R*biaspar)[1]

        if (islog) {
            bias_resp  = rho*(bias_B/B - bias_h/a)
            bias_shift = .
            bias_e     = bias_resp/Ltau
        }
        else {
            bias_shift = Delta*(bias_B/B - bias_h/a)
            bias_resp  = zstar*bias_shift
            bias_e     = bias_shift/(x*Ltau)
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

            bias_resp  = rho*(bias_B/B - bias_h/a)
            bias_shift = .
            bias_e     = bias_resp/Ltau

            return((estimator,bmodel,islog,zstar,t0,t1,tau,lambda,elast,x,rho,Delta, ///
                zlo,zhi,zL,zH,dL,dR,B,bias_h,bias_B,bias_resp,bias_shift,bias_e, ///
                bias_slope,bias_lambda))
        }

        A0 = hi-H
        A1 = (hi^2-H^2)/2
        A2 = (hi^3-H^3)/3

        GtG = J(3,3,0)
        Gtu = J(3,1,0)

        GtG[1..2,1..2] = GtG[1..2,1..2] + pb_intPP(lo,L)

        GtG[1,1] = GtG[1,1] + (1/x^2)*A0
        GtG[1,2] = GtG[1,2] + (1/x^2)*A1
        GtG[2,1] = GtG[1,2]
        GtG[2,2] = GtG[2,2] + (1/x^2)*A2

        GtG[1,3] = GtG[1,3] - (a*A0+m*A1)/x^3
        GtG[3,1] = GtG[1,3]
        GtG[2,3] = GtG[2,3] - (a*A1+m*A2)/x^3
        GtG[3,2] = GtG[2,3]
        GtG[3,3] = GtG[3,3] + (a^2*A0+2*a*m*A1+m^2*A2)/x^4


		gc0 = x*m*(x-1)*zstar
		gc1 = x*m*(x-1)

		Gtu[1] = x*(gc0*A0 + gc1*A1)
		Gtu[2] = x*(gc0*A1 + gc1*A2)
		Gtu[3] = (a*gc0*A0 + (a*gc1 + m*gc0)*A1 + m*gc1*A2)

        biaspar = invsym(GtG)*Gtu

        bias_h = biaspar[1]
        bias_slope = biaspar[2]
        bias_lambda = zstar * (bias_slope/a - (m/a)*(bias_h/a))

        if (bmodel == 1) {
            bias_shift = biaspar[3]
            bias_resp  = zstar*bias_shift
            bias_B     = .
        }
        else {
            R = pb_intP(L,H)/bw
            bias_B = -(R*biaspar[1..2])[1]
            bias_shift = Delta*(bias_B/B - bias_h/a)
            bias_resp  = zstar*bias_shift
        }

        bias_e = bias_shift/(x*Ltau)

        return((estimator,bmodel,islog,zstar,t0,t1,tau,lambda,elast,x,rho,Delta, ///
            zlo,zhi,zL,zH,dL,dR,B,bias_h,bias_B,bias_resp,bias_shift,bias_e, ///
            bias_slope,bias_lambda))
    }

    return(J(1,26,.))
}

end