program polbunchgendata, rclass
	syntax newvarname [, obs(integer 5000) cutoff(real 1) ELasticity(string) ///
		t0(real 0.2) t1(real 0.6) distribution(string) log buncherror(string) ///
		INCOMEeffect(string) PANel(numlist min=3 max=4) ///
		HEAP(numlist min=2 max=3) INATTention(real 0) WELfriction(real 0) ]

	/*
		panel(n T rho [mode]) : generate a pooled panel of n individuals
		observed T years each (obs is overridden to n*T).  All three
		persistence mechanisms leave the marginal distribution of z0 in
		every year EXACTLY distribution() (so A1-A3, E[y_j]=m_j, and the
		closed-form population truth all hold unchanged -- only the
		person-year counts are no longer one multinomial draw):

		  mode 0 (default, or omitted) -- repeat/fresh mixture: each
		    person-year's earnings is a fresh distribution() draw, except
		    that with probability rho it is copied verbatim from the same
		    person's previous year.  corr(z0_it,z0_i,t-k)=rho^k, but the
		    persistence is a two-state jump process (exact repeat, or a
		    fully independent redraw) with no drift in between.  Because the
		    repeats are EXACT, the dependence lands on the diagonal of the
		    bin-count covariance -- vce(conventional) scale(x2) and
		    vce(robust) both recover it.

		  mode 1 (== "smooth") -- Gaussian-copula AR(1): a stationary AR(1)
		    in a latent standard normal (corr(a_it,a_i,t-1)=rho) is mapped
		    through Phi then the triangular inverse-CDF, so a person's
		    earnings drift smoothly through the distribution year to year
		    instead of either freezing or teleporting.  Requires
		    distribution(triangular(a,b,c)).

		  mode 2 (== "components") -- exchangeable Gaussian copula:
		    a_it = sqrt(rho)*xi_i + sqrt(1-rho)*nu_it with xi_i, nu_it iid
		    N(0,1), so every year-pair has the SAME latent correlation rho
		    (a permanent person location xi_i plus a fresh transitory shock
		    each year), mapped through Phi then the triangular inverse-CDF.
		    Unlike mode 0 the person is never in the exact same bin twice,
		    and unlike mode 1 there is no within-person time ordering.  At
		    high rho each person is pinned to a band a few bins wide for all
		    T years: the dependence is spread across the OFF-diagonal of the
		    bin-count covariance, so diagonal-only fixes (scale(x2),
		    vce(robust)) under-state it and only vce(cluster) -- which uses
		    the full co-visitation matrix M -- is calibrated.  Requires
		    distribution(triangular(a,b,c)).

		rho in [0,1): 0 = iid person-years, -> 1 = person frozen across
		years.  Creates the extra variables pid and pyear.
	*/

	/*
		Generates individual-level earnings under the iso-elastic labour
		supply model with a convex kink at cutoff() (rate t0 below, t1
		above).  distribution() gives the COUNTERFACTUAL earnings density
		(what people would earn facing t0 everywhere).

		elasticity(spec)  : the compensated elasticity.  spec is EITHER a
			nonnegative number (homogeneous elasticity) OR any Stata
			expression, evaluated per observation, that draws an individual
			elasticity e_i -- e.g. elasticity(0.4), elasticity(0.4 + rnormal(0,0.1)),
			elasticity(rbeta(2,3)*1.5).  Draws are floored at 1e-8.  The realised
			mean is returned in r(el_mean).  Default 0.4.

		incomeeffect(spec) : income effects.  spec is EITHER a number in
			[0,1) OR a Stata expression, evaluated per observation, drawing
			an individual eta_i in [0,1) (values outside are clipped, with a
			note).  eta is the consumption-curvature parameter of
				u(c,z) = c^(1-eta)/(1-eta) - (n/(1+1/e))(z/n)^(1+1/e).
			elasticity() is the compensated elasticity; the h(z) curvature exponent
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

		inattention(share) / welfriction(theta) : two OPTIMIZATION-FRICTION
			models for would-be bunchers (z0 in the Saez window), both
			nested as special cases of the search-cost model in Chetty,
			Friedman, Olsen & Pistaferri (2011 QJE), "Adjustment Costs,
			Firm Responses, and Micro vs. Macro Labor Supply Elasticities."
			Both are A1 violations: the observed excess mass under-states
			the frictionless behavioural response, so a model that assumes
			frictionless optimization (A1) is misspecified even though A2
			(counterfactual == reference density) and A3 (h0 smooth
			polynomial) still hold. Mutually exclusive; levels only (not
			supported with log or incomeeffect(), which need the
			quasilinear-utility comparison below). z0 here is each
			individual's OWN counterfactual draw (their earnings absent
			the kink) -- Chetty et al.'s "initial offer."

			inattention(share) : Chetty et al.'s Special Case 3 (Section
				II.E) -- a discrete two-type mixture. A fraction `share' of
				the population faces infinite search costs and never
				responds to the kink at all (z stays at z0, exactly as
				drawn); the remaining (1-`share') face none and bunch
				exactly as in the frictionless model. This is the
				mechanism Chetty et al. use to reconcile small "micro"
				bunching elasticities with larger "macro"/structural ones.
				r(n_frictionfail) reports how many would-be bunchers this
				reverted to non-bunchers.

			welfriction(theta) : Chetty et al.'s Special Case 2 (Section
				II.D, eq. 7-8) -- a continuous search-cost/inaction-region
				model, restated as a MONEY-METRIC WELFARE-GAIN threshold
				(theta, a share of income at the cutoff) rather than their
				fixed utils cost kappa, so theta=0.01-0.02 reads as "only
				bunch if it's worth at least 1-2% of income at the kink."
				Under the maintained quasilinear-in-consumption iso-elastic
				disutility v(z)=k/(1+eps)*z^(1+eps), eps=1/e_i, k pinned
				down by z0 solving the flat-t0 FOC, utility differences ARE
				money-metric, so the comparison is exact (no approximation
				needed, unlike Chetty et al.'s quadratic approximation).
				A would-be buncher relocates to the cutoff only if
				u(cutoff) - u(z0) [both evaluated under the TRUE two-rate
				schedule they actually face] exceeds theta*cutoff; workers
				near the edges of the Saez window (where the frictionless
				gain from relocating is smallest) are the first to fail.
				r(n_frictionfail) as above.
	*/

	if "`elasticity'"=="" local elasticity 0.4
	if "`incomeeffect'"=="" local incomeeffect 0

	// constant?  real("...") is nonmissing for a plain number, missing for an expression
	local elc  = real("`elasticity'")
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
			di as error "panel(n T rho [mode]): n and T must both be at least 2."
			exit 198
		}
		if `rhop' < 0 | `rhop' >= 1 {
			di as error "panel(n T rho [mode]): rho must be in [0,1)."
			exit 198
		}
		if !inlist(`smoothp', 0, 1, 2) {
			di as error "panel(n T rho [mode]): mode must be 0 (repeat/fresh mixture), 1 (smooth AR(1) copula) or 2 (exchangeable components copula)."
			exit 198
		}
		if `smoothp' & "`distribution'"!="" & strpos("`distribution'","triangular")!=1 {
			di as error "panel(..., mode `smoothp'): the copula mechanisms (mode 1/2) need distribution(triangular(a,b,c))."
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
			noi di as error "elasticity() must be nonnegative."
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
		if `inattention'<0 | `inattention'>1 {
			noi di as error "inattention() must be in [0,1]."
			exit 198
		}
		if `welfriction'<0 | `welfriction'>=1 {
			noi di as error "welfriction() must be in [0,1)."
			exit 198
		}
		if `inattention'>0 & `welfriction'>0 {
			noi di as error "inattention() and welfriction() are two distinct friction models -- specify only one at a time."
			exit 198
		}
		if (`inattention'>0 | `welfriction'>0) & "`log'"=="log" {
			noi di as error "inattention()/welfriction() are not supported with log (need the levels tax-schedule utility comparison)."
			exit 198
		}
		if (`inattention'>0 | `welfriction'>0) & `iec'>0 {
			noi di as error "inattention()/welfriction() are not supported together with incomeeffect() (need the quasilinear, no-income-effect utility comparison)."
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

		// ---- panel persistence: mixture (0), smooth AR(1) copula (1), exch. copula (2) ----
		if `ispanel' {
			if `smoothp' {
				/*
					Gaussian-copula persistence in a latent standard normal
					a_it (~ N(0,1) marginally every year, so Phi(a_it) is
					U(0,1) and the triangular inverse-CDF below recovers
					distribution() EXACTLY). a,b,c are the triangular params
					parsed above.

					mode 1 -- stationary AR(1):
					  a_it = rho*a_i,t-1 + sqrt(1-rho^2)*e_it,  a_i1 ~ N(0,1)
					  corr(a_it,a_i,t-1) = rho  (smooth year-to-year drift).

					mode 2 -- exchangeable one-factor:
					  a_it = sqrt(rho)*xi_i + sqrt(1-rho)*nu_it
					  corr(a_it,a_i,s) = rho for EVERY t != s  (a permanent
					  person location xi_i + a fresh transitory shock each
					  year; no time ordering, never the exact same value
					  twice). At high rho this pins each person to a narrow
					  band for all T years -> strong positive OFF-diagonal
					  bin-count covariance that scale(x2)/vce(robust) miss.
				*/
				tempvar _lat _u
				sort pid pyear
				if `smoothp' == 1 {
					gen double `_lat' = rnormal() if pyear == 1
					forvalues t = 2/`Tp' {
						by pid: replace `_lat' = `rhop'*`_lat'[_n-1] ///
							+ sqrt(1-`rhop'^2)*rnormal() if pyear == `t'
					}
				}
				else {
					tempvar _xi
					by pid: gen double `_xi' = rnormal() if _n == 1
					by pid: replace `_xi' = `_xi'[1]
					gen double `_lat' = sqrt(`rhop')*`_xi' ///
						+ sqrt(1-`rhop')*rnormal()
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
		cap gen double `ec' = `elasticity'
		if _rc {
			noi di as error "elasticity() must be a nonnegative number or a valid Stata expression, e.g. 0.4 or 0.4 + rnormal(0,0.1)."
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
		// snapshot z0 (each individual's counterfactual draw, pre-response)
		// -- inattention()/welfriction() below need this even after the
		// bunching relocation overwrites `varlist' for bunchers.
		tempvar z0
		gen double `z0' = `varlist'

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

			// ---- optimization frictions (Chetty et al. 2011) ----
			// only would-be bunchers (bunch==1) are at risk of reverting
			// to z0; genuine non-bunchers are untouched either way.
			tempvar frictionfail
			gen byte `frictionfail' = 0
			if `inattention' > 0 {
				// Special Case 3 (Section II.E): discrete two-type
				// mixture -- a fraction `inattention' of the population
				// faces infinite search costs and never responds.
				replace `frictionfail' = (`bunch' & runiform() < `inattention')
			}
			else if `welfriction' > 0 {
				// Special Case 2 (Section II.D, eq. 7-8), restated as a
				// money-metric welfare-gain threshold under quasilinear
				// iso-elastic utility v(z)=k/(1+eps)*z^(1+eps):
				//   - eps_i = 1/e_i
				//   - k_i pinned down by z0_i solving the flat-t0 FOC:
				//       (1-t0) = k_i*z0_i^eps_i  =>  k_i = (1-t0)*z0_i^(-eps_i)
				//   - u(cutoff) uses c(cutoff)=(1-t0)*cutoff (the t0 branch)
				//   - u(z0) uses c(z0)=(1-t0)*cutoff+(1-t1)*(z0-cutoff),
				//     the consumption the individual ACTUALLY gets if they
				//     just keep reporting z0 under the true two-rate
				//     schedule (z0 is in the window, so z0 > cutoff)
				// Relocate to the cutoff only if u(cutoff)-u(z0) clears
				// `welfriction' * cutoff (a share of income at the kink).
				tempvar eps k u_star u_z0 gain
				gen double `eps' = 1/`ec'
				gen double `k'      = (1-`t0') * `z0'^(-`eps') if `bunch'
				gen double `u_star' = (1-`t0')*`cutoff' ///
					- `k'/(1+`eps') * `cutoff'^(1+`eps') if `bunch'
				gen double `u_z0'   = ((1-`t0')*`cutoff' + (1-`t1')*(`z0'-`cutoff')) ///
					- `k'/(1+`eps') * `z0'^(1+`eps') if `bunch'
				gen double `gain'   = `u_star' - `u_z0' if `bunch'
				replace `frictionfail' = (`bunch' & `gain' < `welfriction' * `cutoff')
			}
			replace `varlist' = `z0' if `frictionfail'
			replace `bunch'   = 0    if `frictionfail'

			if `etamax'==0 {
				replace `varlist' = `varlist' * ((1-`t1')/(1-`t0'))^`ec' ///
					if `varlist' > `cutoff' & !`bunch' & !`frictionfail'
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
		local n_frictionfail = 0
		if "`log'" != "log" {
			// `frictionfail' only exists on the levels branch (frictions
			// are not supported with log -- validated above); it is
			// always defined there (initialised to 0) regardless of
			// whether inattention()/welfriction() were actually used.
			count if `frictionfail'
			local n_frictionfail = r(N)
		}
		return scalar n_frictionfail = `n_frictionfail'
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
