*! polbunchbias -- analytical bias of polynomial bunching estimators
*! version 2.0.0  26sep2026
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
        CONstant EXACT SPLITmass POOLmass ]

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

    if "`bw'" == "" local bw = 1
    local _doiter = ("`iterate'" != "")

    // ---- run (optionally iterated) ----
    capture matrix drop _pbxs_core _pbxs_bbeta _pbxs_bcur
    capture scalar drop _pbxs_ecur _pbxs_iter _pbxs_conv _pbxs_diverged
    mata: pbx_bias_solve(`estimator', `bmodel', `islog_val', `useconstant', ///
        `nosplit', `zstar', `t0', `t1', st_matrix("`bcoef'"), `elasticity', ///
        `zlo', `zhi', `zl', `zh', `bw', `_doiter', `tolerance', ///
        `maxiter', `underrelax', "_pbxs")

    tempname out bbeta bcur
    matrix `out'   = _pbxs_core
    matrix `bbeta' = _pbxs_bbeta
    matrix `bcur'  = _pbxs_bcur

    // fitting estimators always carry a fitting-stage bias; a missing
    // bias_h means the degree-`_K' design was singular on the fitted window
    // (estimator 2: or the right window [`zh',`zhi'] alone is too short to
    //  identify the Chetty excess-mass ratio delta)
    if inlist(`estimator',1,2) & missing(el(`out',1,20)) {
        di as error "polbunchbias: the counterfactual polynomial (degree `_K') is not identified on the fitted window [`zlo',`zl'] u [`zh',`zhi'] -- window too short for this order. Lower the polynomial order or widen the estimation window."
    }

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

    di as text _newline "Polynomial bunching bias estimates: estimator `estimator', polynomial order `_K'"
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
    return scalar constant = `useconstant'
    return scalar nosplit = `nosplit'
    return scalar input_elasticity = `elasticity'
    return scalar corrected_elasticity = _pbxs_ecur
    return scalar iterations = _pbxs_iter
    return scalar converged = _pbxs_conv
    return scalar iterate_diverged = _pbxs_diverged

    if `_doiter' & _pbxs_diverged == 1 {
        di as error "Warning: the iterate self-consistency loop did not converge for this cell;"
        di as error "         reporting the un-iterated (non-self-consistent) bias. r(iterate_diverged)=1."
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
// coefficients of p(s + c)   (same degree)
//   new coeff of s^j  =  sum_{k>=j} b_k C(k,j) c^{k-j}
// -------------------------------------------------------------------
real rowvector pbx_polshift(real rowvector b, real scalar c)
{
    real scalar K, j, k
    real rowvector out
    K   = cols(b) - 1
    out = J(1, K+1, 0)
    for (j = 0; j <= K; j++)
        for (k = j; k <= K; k++)
            out[j+1] = out[j+1] + b[k+1]*pbx_choose(k,j)*c^(k-j)
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
// e()-mode: turn the fitted h0 coefficients c (in polbunch's normalised
// coordinate u = (z - zmid)/xscale) into the centred, height-normalised
// polynomial in the centred running variable s that the core wants --
// coefficients of h0hat(s0 + s) where s0 = cest, in s.
//
// Identical for levels and logs: polbunch's `log' means the running
// variable is ALREADY ln(earnings), so `c' is a polynomial in log
// earnings and s = ln z - ln z*.  Exact for any degree K.
// -------------------------------------------------------------------
real rowvector pbx_ecoef_transform(real rowvector c, real scalar cest,
    real scalar xscale)
{
    real rowvector braw
    // h0hat(cest + s/xscale): shift argument by cest, then scale by 1/xscale
    braw = pbx_polscale(pbx_polshift(c, cest), 1/xscale)
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
* The first 26:
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
capture mata: mata drop pbx_bias_solve()

mata:

// reported relative slope: (islog ? 1 : zstar) * h0'(z*)/h0(z*)
real scalar pbx_lambda_report(real rowvector b, real scalar zstar,
                              real scalar islog)
{
    real scalar pref
    if (cols(b) < 2)       return(0)
    if (abs(b[1]) < 1e-14) return(.)
    pref = zstar
    if (islog) pref = 1
    return(pref * b[2] / b[1])
}

real rowvector pbx_bias_core(
    real scalar estimator,
    real scalar bmodel,
    real scalar islog,
    real scalar useconstant,
    real scalar nosplit,
    real scalar zstar,
    real scalar t0,
    real scalar t1,
    real rowvector bcoef,
    real scalar elast,
    real scalar zlo,
    real scalar zhi,
    real scalar zL,
    real scalar zH,
    real scalar bw
)
{
    real scalar tau, Ltau, x, rho, Delta, r, B, K, a0, m0
    real scalar dL, dR, lambda_rep, pref
    real scalar bias_h, bias_B, bias_resp, bias_shift, bias_e
    real scalar bias_slope, bias_lambda
    real scalar atilde, mtilde, Btilde, rtilde, rr, guess
    real scalar lo, hi, L, H, edge_ovh, overhang
    real matrix M, GtG, Glo, Ghi
    real colvector Gtu, Rvec, biaspar, Mhi, dsc, dsc2
    real rowvector ucoef, bhat, h1coef, betat, ures, Rlo, Rhi, Rbar, Jm
    real scalar hminus, hplus, Hstar, Bsaez, sleft, sright
    real scalar sright0, hright0, m_saez, a_saez
    real scalar Asaez, qsaez, disc, xhat, dlogzhat
    real scalar K1, np, delc, q, dqdD, dstep, nlsit, qff
    real scalar uM, Sbar, Sright, lo_right, trueRightMass, trueMass, hlo, llo
    real rowvector biasbeta

    tau   = (1-t0)/(1-t1)
    Ltau  = ln(tau)
    x     = tau^elast
    rho   = ln(x)
    Delta = x - 1
    if (bw <= 0) bw = 1

    K  = cols(bcoef) - 1
    a0 = bcoef[1]
    m0 = 0
    if (K >= 1) m0 = bcoef[2]

    lo = zlo - zstar
    hi = zhi - zstar
    L  = zL  - zstar
    H  = zH  - zstar

    if (islog) r = rho
    else       r = zstar*Delta

    B = pbx_polint(bcoef, 0, r)/bw

    dL = zstar - zL + (zL-zlo)/2
    dR = zH - zstar + (zhi-zH)/2

    lambda_rep = pbx_lambda_report(bcoef, zstar, islog)

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
            return((estimator, bmodel, islog, zstar, t0, t1, tau, lambda_rep,
                elast, x, rho, Delta, zlo, zhi, zL, zH, dL, dR, B,
                J(1, 7, .), J(1, K+1, .)))
        }

        if (islog) ucoef = pbx_polsub(pbx_polshift(bcoef, rho), bcoef)
        else       ucoef = pbx_polsub(pbx_polaffine(bcoef, x, r), bcoef)

        Gtu     = pbx_polmoment(ucoef, H, hi, K)
        biaspar = dsc :* (invsym(M) * (dsc :* Gtu))

        biasbeta = biaspar[|1 \ K+1|]'
        bias_h = biaspar[1]
        if (K >= 1) bias_slope = biaspar[2]

        // reported relative-slope bias -- linearised, exactly as
        // the local-linear polbunchbias.  Pure diagnostic: feeds nothing downstream.
        // (An exact ratio-difference version lands with the .ado rewrite.)
        if (K >= 1 & abs(a0) > 1e-14) {
            pref = zstar
            if (islog) pref = 1
            bias_lambda = pref * (bias_slope/a0 - (m0/a0)*(bias_h/a0))
        }

        // ---- bunching-mass bias ----------------------------------------
        // fit-error piece: -(1/bw) int_E (h0hat - h0) = -(1/bw) int_E P'dbeta
        Rvec   = pbx_polmoment((1), L, H, K)
        bias_B = -(Rvec' * biaspar) / bw

        // overhang / compression term (upper edge only; needs L<0<H)
        if (L < 0 & H > 0) {
            if (islog) edge_ovh = H + rho
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
            bias_e     = bias_resp/Ltau
        }
        else {
            bias_shift = rtilde/zstar - Delta
            bias_e     = ln(1 + rtilde/zstar)/Ltau - elast
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

        if (islog) h1coef = pbx_polshift(betat, rho)
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
            trueMass = pbx_polint(betat, L, H + rho)/bw
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

        // diagonal preconditioner: sc^{-k} on the K+1 beta coefficients,
        // 1 on the (dimensionless) delta index.  Exact change of basis.
        dsc  = pbx_monoscale(lo, hi, L, H, K)
        dsc2 = (dsc \ 1)
        if (pbx_illcond((dsc * dsc') :* (Glo + Ghi))) {
            return((estimator, bmodel, islog, zstar, t0, t1, tau, lambda_rep,
                elast, x, rho, Delta, zlo, zhi, zL, zH, dL, dR, B,
                J(1, 7, .), J(1, K+1, .)))
        }

        delc = Delta
        for (nlsit = 1; nlsit <= 80; nlsit++) {
            q    = 1/(1+delc)
            dqdD = -1/((1+delc)^2)

            // right-region residual: h1_true(s) - q*(P'betat)
            ures = pbx_polsub(h1coef, q :* betat)

            GtG = J(np, np, 0)
            Gtu = J(np, 1, 0)

            GtG[|1,1 \ K1,K1|] = Glo + q^2 :* Ghi
            GtG[|1,np \ K1,np|] = (q*dqdD) :* Mhi
            GtG[|np,1 \ np,K1|] = ((q*dqdD) :* Mhi)'
            GtG[np,np] = dqdD^2 * qff

            Gtu[|1 \ K1|] = q :* pbx_polmoment(ures, H, hi, K)
            Gtu[np] = dqdD * pbx_polint(pbx_polmul(betat, ures), H, hi)

            // missing-mass coefficient on beta: delta/(1+delta) = delc*q,
            // with d/ddelc [delc/(1+delc)] = 1/(1+delc)^2 = q^2
            Jm = J(1, np, 0)
            Jm[|1,1 \ 1,K1|] = Rlo + q :* Rhi + (delc*q) :* Rbar
            Jm[1,np]         = dqdD*Sright + (q^2)*Sbar
            uM = trueMass - ((Rlo + q :* Rhi + (delc*q) :* Rbar) * betat')

            GtG = GtG + Jm' * Jm
            Gtu = Gtu + Jm' * uM

            biaspar = dsc2 :* (invsym((dsc2 * dsc2') :* GtG) * (dsc2 :* Gtu))

            dstep = 0.7 * biaspar[np]
            delc  = delc + dstep
            if (1 + delc <= 1e-6) delc = -1 + 1e-6
            if (abs(dstep) < 1e-11) break
        }

        // The (K+2)-square Gauss-Newton system inverted above is well
        // conditioned at delta = Delta, but its delta row/column collapses
        // (q and dq/ddelta -> 0) once the loop drives delc away -- which it
        // does whenever the right fitting window [zH,zhi] is too short to
        // identify Chetty's delta.  invsym() then silently returns a
        // generalised inverse, the delta step rounds to zero, and the loop
        // "converges" on a garbage delc, yielding a smooth-but-wrong bias
        // branch that kinks into the real one as the window lengthens.  The
        // pre-loop guard cannot see this: it tests only the beta-only Gram
        // Glo + Ghi, which the well-conditioned LEFT window keeps healthy.
        // Re-test the matrix actually inverted and bail (as estimator 1
        // does) when delta is not identified on this window.
        if (pbx_illcond((dsc2 * dsc2') :* GtG)) {
            return((estimator, bmodel, islog, zstar, t0, t1, tau, lambda_rep,
                elast, x, rho, Delta, zlo, zhi, zL, zH, dL, dR, B,
                J(1, 7, .), J(1, K+1, .)))
        }

        biasbeta = biaspar[|1 \ K1|]'
        bias_h = biaspar[1]
        if (K >= 1) bias_slope = biaspar[2]
        if (K >= 1 & abs(a0) > 1e-14) {
            pref = zstar
            if (islog) pref = 1
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
            bias_resp  = zstar*bias_shift
            bias_B     = .
            bias_e     = ln(1 + Delta + bias_shift)/Ltau - elast
        }
        else {
            // reduced-form Bhat = Hstar - int_E h0hat
            Rvec   = pbx_polmoment((1), L, H, K)
            bias_B = -(Rvec' * biaspar[|1 \ K1|]) / bw

            if (L < 0 & H > 0) {
                if (islog) edge_ovh = H + rho
                else       edge_ovh = x*H + r
                overhang = pbx_polint(betat, H, edge_ovh)/bw - B
                if (nosplit) bias_B = bias_B + overhang
                else bias_B = bias_B + overhang ///
                    + (1 - 1/x)*pbx_polint(betat, 0, H)/bw
            }

            bhat   = betat + biaspar[|1 \ K1|]'
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
                bias_e     = bias_resp/Ltau
            }
            else {
                bias_shift = rtilde/zstar - Delta
                bias_e     = ln(1 + rtilde/zstar)/Ltau - elast
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
            if (islog) edge_ovh = H + rho
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
                bias_e     = bias_resp/Ltau
            }
            else {
                bias_shift = rtilde/zstar - Delta
                bias_e     = ln(1 + rtilde/zstar)/Ltau - elast
            }
        }
    }

    if (estimator == 4) {
        // Saez three-region trapezoid.  Intrinsically a TWO-POINT method:
        // it recovers only a line, so its bias against a curved h0 comes
        // from evaluating the exact (curved) h0 over the reference regions.
        // The mass->response transform equation is unchanged -- Saez only
        // ever knows the two reference densities hminus, hplus.
        if (islog) h1coef = pbx_polshift(bcoef, rho)
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
            sright0 = sright + rho
            hright0 = hplus
        }
        else {
            sright0 = (x - 1)*zstar + x*sright
            hright0 = hplus/x
        }
        m_saez = (hright0 - hminus)/(sright0 - sleft)
        a_saez = hminus - m_saez*sleft

        bias_h     = a_saez - a0
        bias_slope = m_saez - m0
        biasbeta = J(1, K + 1, 0)
        biasbeta[1] = bias_h
        if (K >= 1) biasbeta[2] = bias_slope
        pref = zstar
        if (islog) pref = 1
        bias_lambda = pref*(m_saez/a_saez - m0/a0)

        if (islog) {
            if (hminus > 0 & hplus > 0) {
                if (useconstant) dlogzhat = Bsaez*bw/hminus
                else             dlogzhat = 2*Bsaez*bw/(hminus + hplus)
                xhat       = exp(dlogzhat)
                bias_B     = Bsaez - B
                bias_resp  = dlogzhat - rho
                bias_shift = .
                bias_e     = bias_resp/Ltau
            }
        }
        else {
            if (hminus > 0) {
                Asaez = 2*Bsaez*bw/zstar
                if (useconstant) {
                    xhat       = 1 + Asaez/(2*hminus)
                    bias_B     = Bsaez - B
                    bias_resp  = zstar*(xhat - x)
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
                            bias_resp  = zstar*(xhat - x)
                            bias_shift = xhat - x
                            bias_e     = ln(xhat)/Ltau - elast
                        }
                    }
                }
            }
        }
    }

    // columns 1..26 (the local-linear layout); then the per-coefficient bias
    // vector biasbeta (length K+1), so the caller can recover the full
    // fitted polynomial bhat = bcoef + biasbeta.
    return((estimator, bmodel, islog, zstar, t0, t1, tau, lambda_rep, elast,
        x, rho, Delta, zlo, zhi, zL, zH, dL, dR, B, bias_h, bias_B, bias_resp,
        bias_shift, bias_e, bias_slope, bias_lambda, biasbeta))
}

// -------------------------------------------------------------------
// pbx_bias_solve(): optional self-consistent iteration on
// (elasticity, h0 shape), then one final evaluation.  Writes results
// to <stub>_core (1x26), <stub>_bbeta (1xK1), <stub>_bcur (1xK1),
// <stub>_ecur, <stub>_iter, <stub>_conv.
// -------------------------------------------------------------------
void pbx_bias_solve(
    real scalar estimator, real scalar bmodel, real scalar islog,
    real scalar useconstant, real scalar nosplit, real scalar zstar,
    real scalar t0, real scalar t1, real rowvector bhat, real scalar ehat,
    real scalar zlo, real scalar zhi, real scalar zL, real scalar zH,
    real scalar bw, real scalar doiter, real scalar tol, real scalar maxiter,
    real scalar urlx, string scalar stub)
{
    real scalar K1, ecur, enew, iter, conv, dif, k
    real scalar diverged, EBOUND
    real rowvector bcur, bnew, res, bbeta

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

    if (doiter) {
        conv = 0
        for (k = 1; k <= maxiter; k++) {
            res   = pbx_bias_core(estimator, bmodel, islog, useconstant,
                        nosplit, zstar, t0, t1, bcur, ecur,
                        zlo, zhi, zL, zH, bw)
            bbeta = res[|1, 27 \ 1, 26 + K1|]
            if (missing(res[20])) {
                iter = k
                break                        // singular design -- stop
            }
            enew  = ehat - res[24]
            bnew  = bhat - bbeta
            if (abs(bnew[1]) > 1e-14) bnew = bnew :/ bnew[1]   // re-pin height
            // divergence guard: a blown-up coefficient vector or an
            // implausible corrected elasticity means this fixed point is not
            // a contraction here -- stop and fall back to the un-iterated fit.
            if (missing(enew) | hasmissing(bnew) | max(abs(bnew)) > 1e6
                | abs(enew) > EBOUND) {
                iter     = k
                diverged = 1
                break
            }
            dif   = max((abs(enew - ecur), rowmax(abs(bnew - bcur))))
            ecur  = ecur + urlx*(enew - ecur)
            bcur  = bcur + urlx :* (bnew - bcur)
            iter  = k
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

    res   = pbx_bias_core(estimator, bmodel, islog, useconstant, nosplit,
                zstar, t0, t1, bcur, ecur, zlo, zhi, zL, zH, bw)
    if (missing(res[20]) | missing(res[24]) |
        (doiter & abs(ehat - res[24]) > EBOUND)) {
        // the current point gives a missing or runaway bias -- fall back to
        // the raw (fitted) inputs so iterate never degrades the result
        if (doiter & (bcur != bhat | ecur != ehat)) diverged = 1
        bcur = bhat
        ecur = ehat
        conv = 0
        res  = pbx_bias_core(estimator, bmodel, islog, useconstant, nosplit,
                    zstar, t0, t1, bcur, ecur, zlo, zhi, zL, zH, bw)
    }
    bbeta = res[|1, 27 \ 1, 26 + K1|]

    st_matrix(stub + "_core",  res[|1, 1 \ 1, 26|])
    st_matrix(stub + "_bbeta", bbeta)
    st_matrix(stub + "_bcur",  bcur)
    st_numscalar(stub + "_ecur", ecur)
    st_numscalar(stub + "_iter", iter)
    st_numscalar(stub + "_conv", conv)
    st_numscalar(stub + "_diverged", diverged)
}

end
