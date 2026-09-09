program polbunchgendata, rclass
	syntax newvarname [, obs(integer 5000) cutoff(real 1) el(string) ///
		t0(real 0.2) t1(real 0.6) distribution(string) log buncherror(string) ///
		INCOMEeffect(string) PANel(numlist min=3 max=4) ///
		HEAP(numlist min=2 max=3) ]

	/*
		panel(n T rho [smooth]) : generate a pooled panel of n individuals
		observed T years each (obs is overridden to n*T).  Two persistence
		mechanisms, both leaving the marginal distribution of z0 in every
		year EXACTLY distribution() (so A1-A3, E[y_j]=m_j, and the
		closed-form population truth all hold unchanged -- only the
		person-year counts are no longer one multinomial draw):

		  default (smooth=0 or omitted) -- repeat/fresh mixture: each
		    person-year's earnings is a fresh distribution() draw, except
		    that with probability rho it is copied verbatim from the same
		    person's previous year.  corr(z0_it,z0_i,t-k)=rho^k, but the
		    persistence is a two-state jump process (exact repeat, or a
		    fully independent redraw) with no drift in between.

		  smooth=1 -- Gaussian-copula AR(1): a stationary AR(1) in a
		    latent standard normal (corr(a_it,a_i,t-1)=rho) is mapped
		    through Phi then the triangular inverse-CDF, so a person's
		    earnings drift smoothly through the distribution year to year
		    instead of either freezing or teleporting.  Requires
		    distribution(triangular(a,b,c)).

		rho in [0,1): 0 = iid person-years, -> 1 = person frozen across
		years.  Creates the extra variables pid and pyear.
	*/

	/*
		Generates individual-level earnings under the iso-elastic labour
		supply model with a convex kink at cutoff() (rate t0 below, t1
		above).  distribution() gives the COUNTERFACTUAL earnings density
		(what people would earn facing t0 everywhere).

		el(spec)          : the compensated elasticity.  spec is EITHER a
			nonnegative number (homogeneous elasticity) OR any Stata
			expression, evaluated per observation, that draws an individual
			elasticity e_i -- e.g. el(0.4), el(0.4 + rnormal(0,0.1)),
			el(rbeta(2,3)*1.5).  Draws are floored at 1e-8.  The realised
			mean is returned in r(el_mean).  Default 0.4.

		incomeeffect(spec) : income effects.  spec is EITHER a number in
			[0,1) OR a Stata expression, evaluated per observation, drawing
			an individual eta_i in [0,1) (values outside are clipped, with a
			note).  eta is the consumption-curvature parameter of
				u(c,z) = c^(1-eta)/(1-eta) - (n/(1+1/e))(z/n)^(1+1/e).
			el() is the compensated elasticity; the h(z) curvature exponent
			is 1/e_i = 1/ec_i - eta_i.  eta=0 reproduces the
			no-income-effect iso-elastic model exactly.  The bunching window
			is unchanged (width set by the compensated elasticity, as in
			Saez 2010); only the relocation of non-bunchers above the kink
			differs, and it is solved individually.  Not supported with log.
			Default 0.

		buncherror(op expr) : optimisation friction / measurement error for
			bunchers, applied as  cutoff <op expr> , e.g.
			buncherror(+rnormal(0,0.05)) or buncherror(*exp(rnormal(0,0.1))).

		heap(share grid [guardmult]) : round-number heaping.  A random
			`share' of the NON-bunchers report the nearest multiple of
			`grid' instead of their true earnings, EXCEPT within
			guardmult*grid of the cutoff (default guardmult 1.5), which is
			left alone so the excluded window stays clean.  This puts spikes
			in the counterfactual density -- a violation of A3 (h0 a smooth
			polynomial) that leaves A1 and the excess-mass estimand
			E[n_bunchers] intact.  r(n_heaped) reports how many observations
			were moved.  Smaller guardmult puts spikes closer to the window
			(more bias in B); larger keeps B clean but only widens SEs.
	*/

	if "`el'"=="" local el 0.4
	if "`incomeeffect'"=="" local incomeeffect 0

	// constant?  real("...") is nonmissing for a plain number, missing for an expression
	local elc  = real("`el'")
	local iec  = real("`incomeeffect'")

	// -------- panel scaffold --------
	local ispanel = ("`panel'" != "")
	local np  = .
	local Tp  = .
	local rhop = 0
	local smoothp = 0
	if `ispanel' {
		tokenize `panel'
		local np `1'
		local Tp `2'
		local rhop `3'
		if "`4'" != "" local smoothp `4'
		if `np' < 2 | `Tp' < 2 {
			di as error "panel(n T rho [smooth]): n and T must both be at least 2."
			exit 198
		}
		if `rhop' < 0 | `rhop' >= 1 {
			di as error "panel(n T rho [smooth]): rho must be in [0,1)."
			exit 198
		}
		if `smoothp' & "`distribution'"!="" & strpos("`distribution'","triangular")!=1 {
			di as error "panel(..., smooth): the smooth-drift mechanism needs distribution(triangular(a,b,c))."
			exit 198
		}
		local obs = `np' * `Tp'
	}

	quietly {
		clear
		set obs `obs'

		if `ispanel' {
			gen long pid   = ceil(_n / `Tp')
			gen int  pyear = mod(_n - 1, `Tp') + 1
		}

		// -------- validation --------
		if !inrange(`t0',0,1) | !inrange(`t1',0,1) | `t1'<`t0' {
			noi di as error "tax rates t0, t1 must be in [0,1] with t1>t0."
			exit 198
		}
		if `elc'<. & `elc'<0 {
			noi di as error "el() must be nonnegative."
			exit 198
		}
		if `iec'<. & (`iec'<0 | `iec'>=1) {
			noi di as error "incomeeffect() must be in [0,1)."
			exit 198
		}
		if "`log'"=="log" & (`iec'>=. | `iec'!=0) {
			noi di as error "incomeeffect() is not supported with log earnings."
			exit 198
		}

		if "`buncherror'" != "" {
			// validate the expression without touching the (var-less) data
			cap di as text 1*(`cutoff' `buncherror')
			if _rc {
				noi di as error "buncherror() must be an operator + expression, e.g. +rnormal(0,0.05) or *exp(rnormal(0,0.1))."
				exit 198
			}
		}

		// -------- draw counterfactual earnings z0 into `varlist' --------
		tempname r
		if "`distribution'"=="" loc distribution triangular(0,3,0)
		if strpos("`distribution'","triangular")>0 {
			loc _popen  = strpos("`distribution'","(")
			loc _pclose = strpos("`distribution'",")")
			loc _argstr = substr("`distribution'", `_popen'+1, `_pclose'-`_popen'-1)
			loc _argstr = subinstr("`_argstr'", ","," ",.)
			loc _nargs : word count `_argstr'
			if `_nargs' != 3 {
				noi di as error "triangular(a,b,c) requires exactly 3 numeric arguments; got `_nargs' in distribution(`distribution')."
				exit 198
			}
			loc a : word 1 of `_argstr'
			loc b : word 2 of `_argstr'
			loc c : word 3 of `_argstr'
			cap {
				gen double `r'=runiform()
				gen double `varlist'=`a'+sqrt(`r'*(`b'-`a')*(`c'-`a')) if `r'<(`c'-`a')/(`b'-`a')
				replace `varlist'=`b'-sqrt((1-`r')*(`b'-`a')*(`b'-`c')) if `r'>(`c'-`a')/(`b'-`a')
			}
		}
		else cap gen double `varlist'=`distribution'
		if _rc {
			noi di as error "Use a valid Stata random-number function (with parameters) in distribution(); triangular(a,b,c) is also allowed."
			exit 198
		}

		// -------- panel persistence: repeat/fresh mixture, or smooth AR(1) --------
		if `ispanel' {
			if `smoothp' {
				/*
					Gaussian-copula AR(1): stationary latent normal
					a_it = rho*a_i,t-1 + sqrt(1-rho^2)*e_it, e_it~N(0,1),
					a_i1~N(0,1) -- so a_it ~ N(0,1) marginally every year
					and corr(a_it,a_i,t-1)=rho. Phi(a_it) is then U(0,1)
					marginally; pushing it through the triangular inverse-
					CDF recovers distribution() exactly, with SMOOTH drift
					(rather than panel()'s default repeat-or-jump mixture).
					a,b,c are the triangular params parsed above.
				*/
				tempvar _lat _u
				sort pid pyear
				gen double `_lat' = rnormal() if pyear == 1
				forvalues t = 2/`Tp' {
					by pid: replace `_lat' = `rhop'*`_lat'[_n-1] ///
						+ sqrt(1-`rhop'^2)*rnormal() if pyear == `t'
				}
				gen double `_u' = normal(`_lat')
				replace `varlist' = `a' + sqrt(`_u'*(`b'-`a')*(`c'-`a')) ///
					if `_u' < (`c'-`a')/(`b'-`a')
				replace `varlist' = `b' - sqrt((1-`_u')*(`b'-`a')*(`b'-`c')) ///
					if `_u' >= (`c'-`a')/(`b'-`a')
			}
			else if `rhop' > 0 {
				sort pid pyear
				forvalues t = 2/`Tp' {
					by pid: replace `varlist' = `varlist'[_n-1] ///
						if pyear == `t' & runiform() < `rhop'
				}
			}
		}

		su `varlist', meanonly
		if !inrange(`cutoff', r(min), r(max)) {
			noi di as error "cutoff=`cutoff' is not within the support of distribution()."
			exit 301
		}

		// -------- individual compensated elasticity --------
		tempvar ec
		cap gen double `ec' = `el'
		if _rc {
			noi di as error "el() must be a nonnegative number or a valid Stata expression, e.g. 0.4 or 0.4 + rnormal(0,0.1)."
			exit 198
		}
		count if `ec'<0 | missing(`ec')
		if r(N) {
			noi di as text "note: `r(N)' of `obs' elasticity draws were negative or missing; floored at 1e-8."
		}
		replace `ec' = max(`ec', 1e-8)

		// -------- individual income-effect parameter eta --------
		tempvar etav
		cap gen double `etav' = `incomeeffect'
		if _rc {
			noi di as error "incomeeffect() must be a number in [0,1) or a valid Stata expression."
			exit 198
		}
		count if `etav'<0 | `etav'>=1 | missing(`etav')
		if r(N) {
			noi di as text "note: `r(N)' of `obs' incomeeffect draws were outside [0,1) or missing; clipped."
		}
		replace `etav' = max(0, min(`etav', 1-1e-6))
		su `etav', meanonly
		local etamax = r(max)

		// -------- apply the behavioural response --------
		tempvar bunch

		if "`log'"=="log" {
			tempvar rho
			gen double `rho' = `ec' * (ln(1-`t0') - ln(1-`t1'))
			gen byte `bunch' = inrange(`varlist', `cutoff', `cutoff' + `rho')
			replace `varlist' = `cutoff' if `bunch'
			replace `varlist' = `varlist' - `rho' if `varlist' > `cutoff' & !`bunch'
		}
		else {
			// bunching window: z0 in (cutoff, cutoff * ((1-t0)/(1-t1))^ec]
			// (this holds with or without income effects -- Saez 2010)
			gen byte `bunch' = inrange(`varlist', `cutoff', ///
				`cutoff' * ((1-`t0')/(1-`t1'))^`ec')
			replace `varlist' = `cutoff' if `bunch'

			if `etamax'==0 {
				replace `varlist' = `varlist' * ((1-`t1')/(1-`t0'))^`ec' ///
					if `varlist' > `cutoff' & !`bunch'
			}
			else {
				mata: pbgd_income("`varlist'", "`ec'", "`etav'", "`bunch'", ///
					`t0', `t1', `cutoff')
			}
		}

		// -------- optimisation friction / measurement error for bunchers --------
		if "`buncherror'" != "" {
			replace `varlist' = `cutoff' `buncherror' if `bunch'
		}

		// -------- round-number heaping among NON-bunchers -----------------
		//  heap(share grid): a random `share' of non-bunchers report the
		//  nearest multiple of `grid'.  Deterministic mean shift (spikes in
		//  h0), not extra noise -- an A3 violation (h0 not a smooth
		//  polynomial) that leaves A1 and the excess-mass estimand intact.
		//  Heaping is suppressed within 1.5*grid of the cutoff so the
		//  excluded region and the behavioural window stay clean; choose
		//  `grid' so the heap points fall in the reference region.
		local heap_n = 0
		if "`heap'" != "" {
			gettoken _hshare _hrest : heap
			gettoken _hgrid  _hguard : _hrest
			if "`_hguard'" == "" local _hguard 1.5
			if `_hgrid' <= 0 {
				noi di as error "heap(share grid [guardmult]): grid must be positive."
				exit 198
			}
			if `_hshare' < 0 | `_hshare' > 1 {
				noi di as error "heap(share grid [guardmult]): share must be in [0,1]."
				exit 198
			}
			tempvar _hpick
			gen byte `_hpick' = !`bunch' & runiform() < `_hshare' ///
				& abs(`varlist' - `cutoff') >= `_hguard'*`_hgrid'
			replace `varlist' = round(`varlist'/`_hgrid')*`_hgrid' if `_hpick'
			count if `_hpick'
			local heap_n = r(N)
		}

		count if `bunch'
		return scalar n_bunchers     = r(N)
		return scalar share_bunching = r(N)/`obs'
		return scalar n_heaped       = `heap_n'
		su `ec', meanonly
		return scalar el_mean        = r(mean)
		su `etav', meanonly
		return scalar incomeeffect   = r(mean)
		return scalar obs            = `obs'
		if `ispanel' {
			return scalar n_panel      = `np'
			return scalar T_panel      = `Tp'
			return scalar rho_panel    = `rhop'
			return scalar smooth_panel = `smoothp'
		}
	}
end

// ---------------------------------------------------------------------
//  Income-effects relocation.  For each non-buncher above the kink,
//  solve the upper-segment first-order condition for the relocation
//  ratio zeta = z1/z0.  Eliminating ability n, the FOC reduces to
//
//     tau * ( tau*zeta + (1-tau)*w )^(-eta_i)  =  zeta^(1/e_i),
//
//  with tau = (1-t1)/(1-t0),  w = zstar/z0,  1/e_i = 1/ec_i - eta_i.
//  G(zeta) is strictly decreasing; the root lies in (w, 1) for a
//  non-buncher (G(w) >= 0, G(1) < 0).  eta_i = 0 gives zeta = tau^ec_i,
//  the iso-elastic proportional relocation.  eta_i varies across
//  individuals (heterogeneous income effects).
// ---------------------------------------------------------------------
capture mata: mata drop pbgd_income()
mata:
void pbgd_income(string scalar zv, string scalar ecv, string scalar etav,
                 string scalar bv, real scalar t0, real scalar t1,
                 real scalar zstar)
{
	real scalar tau, n, i, w, ie, eta, zlo, zhi, zm, g, it
	real colvector z, ec, et, b

	st_view(z=., ., zv)
	ec = st_data(., ecv)
	et = st_data(., etav)
	b  = st_data(., bv)
	tau = (1-t1)/(1-t0)
	n   = rows(z)

	for (i=1; i<=n; i++) {
		if (b[i] | z[i] <= zstar) continue

		eta = et[i]
		ie  = 1/ec[i] - eta          // 1/e_i, the h(z) curvature exponent
		if (ie < 1e-4) ie = 1e-4
		w   = zstar / z[i]           // in (0,1)

		zlo = w ; zhi = 1
		for (it=1; it<=200; it++) {
			zm = 0.5*(zlo+zhi)
			g  = tau*(tau*zm + (1-tau)*w)^(-eta) - zm^ie
			if (g > 0) zlo = zm
			else       zhi = zm
			if (zhi-zlo < 1e-13) break
		}
		z[i] = 0.5*(zlo+zhi) * z[i]
	}
}
end
