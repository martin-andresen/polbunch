program bunchcalc, rclass
    syntax [anything], estimator(integer) cutoff(real) polynomial(integer) bw(real) ///
        z(name) h(integer) l(integer) ///
        [t0(numlist min=1 max=1 <1) t1(numlist min=1 max=1 <1) ///
         boot log nonormalize constant]

    local nlcom `"`anything'"'

    local cutoff_orig = `cutoff'
    local bw_orig     = `bw'

    if "`nonormalize'" == "nonormalize" {
        local cutoff_est = `cutoff_orig'
        local bw_est     = `bw_orig'
        local xscale     = 1
    }
    else {
        local cutoff_est = 0
        local bw_est     = 1
        local xscale     = `bw_orig'
    }

    // polynomial term names
    forvalues k = 1/`polynomial' {
        local term`k' c.`z'
        if `k' > 1 {
            forvalues j = 2/`k' {
                local term`k' `term`k''#c.`z'
            }
        }
    }

    // h0(cutoff)
    local h0cut _b[h0:_cons]
    forvalues k = 1/`polynomial' {
        local h0cut `h0cut' + _b[h0:`term`k'']*`cutoff_est'^`k'
    }

    // B expression
    if inlist(`estimator', 0, 1) {
        local B 0
        forvalues j = 1/`=`h'+`l'' {
            local B `B' + _b[bunching:`j'.bunch]
        }
    }
    else {
        local ex_lo = `cutoff_est' + (-`l' + 1)*`bw_est'
        local ex_hi = `cutoff_est' + (`h' + 1)*`bw_est'

        local cfmass _b[h0:_cons]*(`ex_hi' - `ex_lo')
        forvalues k = 1/`polynomial' {
            local cfmass `cfmass' + ///
                _b[h0:`term`k'']*((`ex_hi'^(`k'+1) - `ex_lo'^(`k'+1))/(`k'+1))
        }

        local B (_b[bunching:Hstar] - (`cfmass'))
    }

    // Build nlcom only if not supplied
    if `"`nlcom'"' == "" {
        local nlcom

        // h0 coefficients
        forvalues k = 1/`polynomial' {
            local nlcom `nlcom' (h0_`k': _b[h0:`term`k''])
        }
        local nlcom `nlcom' (h0_cons: _b[h0:_cons])

        // h1 coefficients
        forvalues k = 1/`polynomial' {
            if `estimator' == 0 {
                local h1expr _b[h1:`term`k'']
            }
            else if `estimator' == 1 {
                local h1expr _b[h0:`term`k'']
            }
            else if `estimator' == 2 {
                local h1expr (_b[h0:`term`k'']/(1+_b[bunching:delta]))
            }
            else if "`log'" == "" {
                local h1expr (_b[h0:`term`k'']*(1+_b[bunching:delta])^(`k'+1))
            }
            else {
                local h1expr _b[h0:`term`k'']
                if `polynomial' > `k' {
                    forvalues n = `=`k'+1'/`polynomial' {
                        local h1expr (`h1expr' + ///
                            _b[h0:`term`n'']*comb(`n',`k')*ln(1+_b[bunching:delta])^(`n'-`k'))
                    }
                }
            }
            local nlcom `nlcom' (h1_`k': `h1expr')
        }

        // h1 constant
        if `estimator' == 0 {
            local h1cons _b[h1:_cons]
        }
        else if `estimator' == 1 {
            local h1cons _b[h0:_cons]
        }
        else if `estimator' == 2 {
            local h1cons (_b[h0:_cons]/(1+_b[bunching:delta]))
        }
        else if "`log'" == "" {
            local h1cons (_b[h0:_cons]*(1+_b[bunching:delta]))
        }
        else {
            local h1cons _b[h0:_cons]
            forvalues n = 1/`polynomial' {
                local h1cons (`h1cons' + _b[h0:`term`n'']*ln(1+_b[bunching:delta])^`n')
            }
        }
        local nlcom `nlcom' (h1_cons: `h1cons')

        // common bunching objects
        local nlcom `nlcom' ///
            (number_bunchers: `B') ///
            (excess_mass: (`B')/(`h0cut'))

        // shift/MR/elasticity: estimator 3 directly use delta
        if inlist(`estimator', 3) {
            local nlcom `nlcom' ///
                (shift: _b[bunching:delta]) ///
                (marginal_response: _b[bunching:delta]*`cutoff_orig')

            if "`t0'" != "" & "`t1'" != "" {
                local nlcom `nlcom' ///
                    (elasticity: ln(1+_b[bunching:delta])/(ln(1-(`t0'))-ln(1-(`t1'))))
            }
        }

        // constant approximation for estimator 0/1/2: algebraic, so use nlcom
        if inlist(`estimator', 0, 1,2) & "`constant'" != "" {
            local deltaz ((`B')*`bw_orig'/(`h0cut'))

            if "`log'" == "" {
                local shift_expr (`deltaz'/`cutoff_orig')
                local mr_expr    (`deltaz')
                local el_expr    (ln(1+`deltaz'/`cutoff_orig')/(ln(1-(`t0'))-ln(1-(`t1'))))
            }
            else {
                local shift_expr (exp(`deltaz')-1)
                local mr_expr    (`deltaz')
                local el_expr    (`deltaz'/(ln(1-(`t0'))-ln(1-(`t1'))))
            }

            local nlcom `nlcom' ///
                (shift: `shift_expr') ///
                (marginal_response: `mr_expr')

            if "`t0'" != "" & "`t1'" != "" {
                local nlcom `nlcom' (elasticity: `el_expr')
            }
        }
    }

    tempname b V

    // Evaluate algebraic transformations
    if "`boot'" == "" {
        quietly nlcom `nlcom'
        matrix `b' = r(b)
        matrix `V' = r(V)
    }
    else {
        local nlcom2 `"`nlcom'"'
        while `"`nlcom2'"' != "" {
            gettoken use nlcom2 : nlcom2, match(parns)
            matrix `b' = nullmat(`b'), `=`use''
        }
    }

    // Estimator 0/1/2, non-constant: append eresp-based shift/MR/elasticity
    if inlist(`estimator', 0, 1,2) & "`constant'" == "" {
        tempname h0coef

        matrix `h0coef' = J(1, `=`polynomial'+1', .)
        forvalues k = 1/`polynomial' {
            matrix `h0coef'[1,`k'] = _b[h0:`term`k'']
        }
        matrix `h0coef'[1,`=`polynomial'+1'] = _b[h0:_cons]

        mata: st_numscalar("rhat", eresp(`=`B'', `cutoff_est', st_matrix("`h0coef'"), `bw_est', `xscale'))

        if rhat == . {
            return scalar exit = 1
            exit
        }

        if "`log'" == "" {
            matrix `b' = `b', `=rhat/`cutoff_orig'', rhat
        }
        else {
            matrix `b' = `b', `=exp(rhat)-1', rhat
        }

        if "`t0'" != "" & "`t1'" != "" {
            if "`log'" == "" {
                matrix `b' = `b', `=ln(1+rhat/`cutoff_orig')/(ln(1-(`t0'))-ln(1-(`t1')))'
            }
            else {
                matrix `b' = `b', `=rhat/(ln(1-(`t0'))-ln(1-(`t1')))'
            }
        }

        // Manual delta-method VCV for appended objects
        if "`boot'" == "" {
            tempname Graw Gextra Vraw Vextra
            matrix `Vraw' = e(V)
            matrix `Graw' = J(1, colsof(e(b)), 0)

            local r_est = rhat/`xscale'
            local zu = `cutoff_est' + `r_est'

            local hzu _b[h0:_cons]
            forvalues k = 1/`polynomial' {
                local hzu `hzu' + _b[h0:`term`k'']*(`zu'^`k')
            }

            // dr/dB in original units
            local drdB = `xscale'*`bw_est'/(`hzu')

            // dr/dbeta_k in original units
            forvalues k = 1/`polynomial' {
                local intk = ((`zu'^(`k'+1) - `cutoff_est'^(`k'+1))/(`k'+1))
                local col = colnumb(e(b), "h0:`term`k''")
                matrix `Graw'[1,`col'] = -`xscale'*`intk'/(`hzu')
            }

            local int0 = (`zu' - `cutoff_est')
            local col = colnumb(e(b), "h0:_cons")
            matrix `Graw'[1,`col'] = -`xscale'*`int0'/(`hzu')

            // B=sum bunch dummies
            forvalues j = 1/`=`l'+`h'' {
                local col = colnumb(e(b), "bunching:`j'.bunch")
                matrix `Graw'[1,`col'] = `drdB'
            }

            local nextra = 2
            if "`t0'" != "" & "`t1'" != "" local nextra = 3

            matrix `Gextra' = J(`nextra', colsof(e(b)), 0)

            if "`log'" == "" {
                matrix `Gextra'[1,1] = `Graw'/`cutoff_orig'
                matrix `Gextra'[2,1] = `Graw'
                if `nextra' == 3 {
                    local A = ln(1-(`t0')) - ln(1-(`t1'))
                    matrix `Gextra'[3,1] = (1/((`cutoff_orig' + rhat)*`A')) * `Graw'
                }
            }
            else {
                matrix `Gextra'[1,1] = exp(rhat)*`Graw'
                matrix `Gextra'[2,1] = `Graw'
                if `nextra' == 3 {
                    local A = ln(1-(`t0')) - ln(1-(`t1'))
                    matrix `Gextra'[3,1] = (1/`A') * `Graw'
                }
            }

            matrix `Vextra' = `Gextra' * `Vraw' * `Gextra''

            matrix `V' = ///
                (`V', J(rowsof(`V'), `nextra', 0) \ ///
                 J(`nextra', colsof(`V'), 0), `Vextra')
        }
    }

    // Names
    local cnames

    forvalues k = 1/`polynomial' {
        local cnames `cnames' h0:`term`k''
    }
    local cnames `cnames' h0:_cons

    forvalues k = 1/`polynomial' {
        local cnames `cnames' h1:`term`k''
    }
    local cnames `cnames' h1:_cons

    local cnames `cnames' ///
        bunching:number_bunchers ///
        bunching:excess_mass
	
	if `estimator'==2 {
		local cnames `cnames' bunching:delta
		}
		
	 local cnames `cnames' ///
        bunching:shift ///
        bunching:marginal_response

    if "`t0'" != "" & "`t1'" != "" {
        local cnames `cnames' bunching:elasticity
    }

    matrix colnames `b' = `cnames'
    if "`boot'" == "" {
        matrix colnames `V' = `cnames'
        matrix rownames `V' = `cnames'
    }

    return local nlcom `"`nlcom'"'
    return matrix b = `b'
    if "`boot'" == "" return matrix V = `V'
    return scalar exit = 0
end

mata:

real scalar eresp(real scalar B, real scalar tau, real matrix cf, real scalar bw, real scalar xscale)
{
    real matrix integral, roots, realroots, out

    if (cols(cf) == 1) {
        if (cf[1] <= 0) return(.)
        return(((B * bw) / cf[1]) * xscale)
    }

    integral = polyinteg(cf, 1)
    integral[1] = -polyeval(integral, tau) - B*bw

    roots = polyroots(integral)
    realroots = Re(select(roots, Im(roots) :== 0))
    out = sort(select(realroots, realroots :> tau)', 1)'

    if (cols(out) == 0) return(.)
    else return((out[1] - tau) * xscale)
}
end

