				*! polbunch version date 20260910
				* Author: Martin Eckhoff Andresen
				* This program is part of the polbunch package.
				
				cap prog drop polbunch
				program polbunch, eclass sortpreserve
					syntax varlist(min=1 max=2) [if] [in],  CUToff(real) [bw(numlist min=1 max=1 >0)  ///
					LIMits(numlist max=2 min=2 integer) ///
					t0(numlist min=1 max=1 <1) ///
					t1(numlist min=1 max=1 <1) ///
					POLynomial(integer 7) ///
					NOIsily ///
					ESTimator(integer 3) /// Specify estimator - 3 = theoretically consistent efficient estimator, 2 = chetty, 1 = no adjustment, 0=data to the left only, 4=Saez three-region trapezoid approximaton
					DELTAmax(real 1) /// upper bound on the structural delta in the estimator 2/3 profile search (default 1 = a 100% earnings response)
					ALLOWnegative /// allow the structural shift delta to range down to -1 in the estimator 2/3 profile; by default the search is restricted to delta >= 0 (a convex kink implies a non-negative response)
					nodrop ///
					notransform ///
					positive ///
					nonormalize ///
					BOOTreps(integer 500) ///
					vce(string) ///
					SCALE(string) ///
					NOMASScorr ///
					log ///
					constant ///
					EXACT ///
					POOLmass ///
					SPLITmass ///
					nodots /// suppress dots for bootstrap progress
					test(string) ///
					nozero ///
					norankred ///
					noRANKcheck ///
					NOBias ///
					NOITERate ///
					SAVEbins(string) ///
					CONTRAST ///
					]

					local cmdline0 `"`0'"'

					 quietly {
							if "`t0'"!="" {
							if "`t1'"=="" {
								noi di as error "If specifying one tax rate (options t0 or t1), specify both."
								exit 301
							}
							if `t0'==`t1' {
								noi di as error "Tax rates t0 and t1 cannot be equal - no incentive! Estimate reduced form bunching by omitting tax rates."
								exit 301
							}
							if `t1'<`t0' {
								noi di as error "Polbunch currently support only convex kinks, t1>t0."
								exit 301
							}
						}
						
						
						if "`vce'"!="none" {
							loc coeftabresults=c(coeftabresults)
							set coeftabresults off
						}
						if !inlist(`estimator',0,1,2,3,4) {
							noi di as error "Option estimator can take only values 0 (using data to the left only),  1 (no adjustment), 2 (Chetty et. al. adjustment),  3 (theoretically consistent and efficient estimator) or 4 (Saez trapezoid approximation)."
							exit 301
						}

						/* varlist width: 1 var (+bw) = individual data, 2 vars = pre-binned.
						   Needed here already so vce(cluster ...) can tell whether its
						   token is a cluster variable (individual) or a co-visitation
						   stub (pre-binned). */
						loc nvars : word count `varlist'

						/*
							Two toggle axes, mirroring polbunchbias:

							  constant / exact -- how excess mass is turned into a
							    response (constant-density approximation vs. exact
							    inversion of the counterfactual-density integral).
							    Default: constant for estimators 1 and 2, exact for
							    estimators 0 and 3.

							  poolmass / splitmass -- how the bunching mass B is
							    formed from the observed excluded-region mass M:
							      poolmass  : B = M - int_{zL}^{zH} h0
							      splitmass : B = M - int_{zL}^{z*} h0
							                       - int_{z*}^{zH} h1
							    Default: poolmass for estimators 1 and 2, splitmass
							    for estimators 0 and 3.  Irrelevant for estimator 1
							    (h0 == h1) and ignored for estimator 4 (Saez has its
							    own two-point mass calculation).
						*/
						if "`constant'" != "" & "`exact'" != "" {
							noi di as error "Specify at most one of constant / exact."
							exit 198
						}
						if "`poolmass'" != "" & "`splitmass'" != "" {
							noi di as error "Specify at most one of poolmass / splitmass."
							exit 198
						}

						/*
							Sign of the structural shift.  polbunch supports only
							convex kinks (t1 > t0), where the marginal buncher
							reduces earnings to the cutoff -- so delta >= 0 by
							construction.  The estimator 2/3 profile is therefore
							restricted to delta >= 0 by default (equivalent to the
							old `positive' option, which is now redundant).
							`allownegative' opens the search back up to delta > -1
							for the rare finite-sample case where the true delta
							is ~0 and noise pulls the profile minimum slightly
							below zero.
						*/
						if "`positive'" != "" & "`allownegative'" != "" {
							noi di as error "Specify at most one of positive / allownegative."
							exit 198
						}
						if "`allownegative'" == "" local positive positive

						if "`constant'" != ""      loc useconstant = 1
						else if "`exact'" != ""    loc useconstant = 0
						else                       loc useconstant = inlist(`estimator',1,2)

						if "`poolmass'" != ""      loc nosplit = 1
						else if "`splitmass'" != "" loc nosplit = 0
						else                       loc nosplit = inlist(`estimator',1,2)
						if `estimator' == 1        loc nosplit = 1

						loc constant = cond(`useconstant',"constant","")

						if "`test'"=="" {
							if `estimator'==4 loc test wald
							else              loc test all
						}
						else {
							if !inlist("`test'","none","minimumdistance","wald","hausman","all","forceall") {
								noi di as error "test() can only take wald, minimumdistance, hausman, all or none."
								exit 301
							}
							if `estimator'==4 & !inlist("`test'","wald","all","none") {
								noi di as error "Only test(wald) supported for estimator 4 (Saez) - simple linear restrictions."
								exit 301
							}
							if `estimator'==1 & "`test'"=="minimumdistance" {
								noi di as error "For the naive estimator test(wald) IS the minimum-distance test (no nuisance parameter). Use test(wald) or test(hausman)."
								exit 301
							}
							if "`test'"=="wald" & !inlist(`estimator',1,4) {
								noi di as txt "note: Wald test is a conditional shape diagnostic at the mass-implied delta;"
								noi di as txt "      minimumdistance or hausman is recommended for formal testing."
							}
						}
						
						

						if `estimator'==4 loc polynomial=0

						tempvar touse
						marksample touse
						preserve
						drop if !`touse'
						
						/*
							vce() -- Stata-standard form  vce(type [, subopts]).

							  type
							    conventional  (default; synonyms analytic /
							        unadjusted / oim) -- multinomial-count delta
							        method, Var(y_j) ~ multinomial ~ Poisson.
							    robust / hc0 / hc1 / hc2 / hc3 -- Eicker-White
							        residual meat diag((y_j - yhat_j)^2) on the
							        polynomial-fit rows; robust = hc1.  (hc0/hc1
							        accepted but undocumented -- regress exposes
							        only robust/hc2/hc3.)
							    bootstrap [, multinomial | residual | wild
							                 bayesian  normal | bc | percentile
							                 wildweights() reps() seed() ]
							        multinomial (default): Dirichlet resample of the
							        bin counts (= individual resample).  residual:
							        Chetty/CFOP residual bootstrap over bins.  wild:
							        heteroskedasticity-robust wild bootstrap over the
							        bin residuals (resampling twin of vce(robust)).
							    none -- no inference.

							Internal legacy locals kept:
							  `vce'      = analytic | bootstrap | none (control flow)
							  hctype     = -1 multinomial ; 0..3 HC flavour
							  boottype   = multinomial | bayesian | residual | wild
							  bootci     = normal | bc | percentile
							  wildwt     = rademacher | mammen | webb
							  scalemode  = 1 | x2 | <#>   (glm-style, conventional only)
							  masscorr   = 1 unless nomasscorr
						*/
						gettoken vtype vsub : vce, parse(",")
						local vtype = strtrim("`vtype'")
						/* pull the first whitespace-delimited token; the remainder
						   is the cluster variable for vce(cluster clustvar) */
						gettoken vtype _vclhead : vtype, parse(" ")
						local vtype   = strlower(strtrim("`vtype'"))
						local clustvar = strtrim("`_vclhead'")
						gettoken _vc vsub : vsub, parse(",")
						local vsub = strtrim(`"`vsub'"')

						if "`vtype'"=="" | inlist("`vtype'","analytic","unadjusted","oim") ///
							loc vtype conventional
						if "`vtype'"=="bayes" {
							loc vtype bootstrap
							loc vsub `vsub' bayesian
						}

						if !inlist("`vtype'","conventional","robust","hc0","hc1","hc2") ///
							& !inlist("`vtype'","hc3","bootstrap","none","cluster") {
							noi di as error "vce() type must be conventional, robust, hc2, hc3, cluster, bootstrap or none."
							exit 198
						}

						local nclusters .
						if "`vtype'"=="cluster" {
							if `nvars'==1 {
								/* individual data: the token is the cluster variable;
								   M and the cluster count are built internally. */
								if "`clustvar'"=="" {
									noi di as error "vce(cluster) requires a cluster variable, e.g. vce(cluster id)."
									exit 198
								}
								capture confirm variable `clustvar'
								if _rc {
									noi di as error "vce(cluster `clustvar'): variable `clustvar' not found."
									exit 111
								}
								if `"`vsub'"' != "" {
									noi di as error "vce(cluster `clustvar') takes no suboptions with individual-level data (nclusters() is for pre-binned co-visitation input)."
									exit 198
								}
								quietly count if missing(`clustvar')
								if r(N) > 0 quietly drop if missing(`clustvar')
							}
							else {
								/* pre-binned data: the token is the co-visitation
								   stub -- either J variables stub1..stubJ (full
								   J x J matrix, row-aligned with the histogram) or a
								   single variable stub (diagonal-only form).  The
								   number of clusters must be given explicitly. */
								if "`clustvar'"=="" {
									noi di as error "vce(cluster) on pre-binned data needs a co-visitation stub, e.g. vce(cluster m, nclusters(#))."
									exit 198
								}
								_pb_parse_clustersub , `vsub'
								local nclusters `s(nclusters)'
								if "`nclusters'"=="." {
									local nclusters : char _dta[polbunch_nclusters]
									if "`nclusters'"=="" local nclusters .
								}
								if "`nclusters'"=="." {
									noi di as error "vce(cluster `clustvar') on pre-binned data requires nclusters(#), the number of clusters (individuals)."
									exit 198
								}
								if `nclusters' < 2 {
									noi di as error "nclusters() must be at least 2."
									exit 198
								}
							}
						}
						else local clustvar ""

						loc hctype = -1
						if inlist("`vtype'","robust","hc1") loc hctype = 1
						else if "`vtype'"=="hc0" loc hctype = 0
						else if "`vtype'"=="hc2" loc hctype = 2
						else if "`vtype'"=="hc3" loc hctype = 3

						loc boottype multinomial
						loc bootci   normal
						loc wildwt   rademacher
						if "`vtype'"=="bootstrap" & `"`vsub'"'!="" {
							/*
								Parse the vce(bootstrap, <subopts>) list in a
								sub-program so its own -syntax- does not clobber
								the caller's varlist / main options.
							*/
							_pb_parse_vcesub , `vsub'
							loc boottype `s(boottype)'
							loc bootci   `s(bootci)'
							if `"`s(wildwt)'"'!="" {
								loc wildwt = strlower("`s(wildwt)'")
								if !inlist("`wildwt'","rademacher","mammen","webb") {
									noi di as error "wildweights() must be rademacher, mammen or webb."
									exit 198
								}
							}
							if `s(reps)'>1 loc bootreps `s(reps)'
							if `"`s(seed)'"'!="" set seed `s(seed)'
						}

						if "`vtype'"=="none"           loc vce none
						else if "`vtype'"=="bootstrap"  loc vce bootstrap
						else                            loc vce analytic

						local grad = cond("`vce'"=="analytic","","nograd")

						/*
							contrast -- paired bootstrap test of this estimator's
							elasticity against the model-consistent efficient
							reference (estimator 3, exact inversion, splitmass),
							refit on the SAME resampled histogram every replication
							so the two estimates' dependence is carried directly.
							Same reference bins and window for both.  Reported in
							the model-test block; stored in e(contrast_*).
						*/
						local docontrast = ("`contrast'" != "")
						if `docontrast' {
							if "`vce'" != "bootstrap" {
								noi di as error "contrast requires vce(bootstrap[, ...]) -- the test is a paired bootstrap of the elasticity difference."
								exit 198
							}
							if !inlist(`estimator',1,2,3) {
								noi di as error "contrast is available for estimator(1), (2) or (3) only -- estimators 0 and 4 use a different reference-bin set / window than the estimator(3) reference."
								exit 198
							}
							if "`t0'"=="" | "`t1'"=="" {
								noi di as error "contrast tests the elasticity, so tax rates t0() and t1() are required."
								exit 198
							}
							if "`transform'"=="notransform" {
								noi di as error "contrast is incompatible with notransform (it needs the elasticity)."
								exit 198
							}
							if `estimator'==3 & `useconstant'==0 & `nosplit'==0 {
								noi di as error "the reported estimate already IS the reference (estimator(3), exact, splitmass) -- nothing to contrast."
								noi di as error "Use estimator(3) with constant or poolmass, or estimator(1)/(2), to contrast against it."
								exit 198
							}
						}

						/* scale() -- glm-style overdispersion multiplier; conventional only */
						if "`scale'"=="" loc scalemode 1
						else {
							loc scalemode = strlower(strtrim("`scale'"))
							if !inlist("`scalemode'","1","x2") & missing(real("`scalemode'")) {
								noi di as error "scale() must be 1, x2 or a positive number."
								exit 198
							}
							if !missing(real("`scalemode'")) {
								if real("`scalemode'")<=0 {
									noi di as error "scale() must be positive."
									exit 198
								}
							}
							if "`scalemode'"!="1" & "`vtype'"!="conventional" {
								noi di as error "scale() is supported only with vce(conventional)."
								exit 198
							}
						}
						loc masscorr = ("`nomasscorr'"=="")
						
						if `bootreps'<=1 {
							noi di as error "Option bootreps can only take an integer >1 (binned bootstrap)."
							exit 301
						}
						if `polynomial'<0 {
							noi di as error "Polynomial must be a nonnegative integer"
							exit 301
						}
						if `deltamax'<=0 {
							noi di as error "Option deltamax() must be a positive number."
							exit 301
						}
						
						
						// check varlist vs bw opts
						if (`nvars'==1&"`bw'"=="")|(`nvars'==2&"`bw'"!="") {
							noi di as error "Varlist must either contain 1 variable (earnings z) and option bw be specified (individual level data) or 2 variables (frequency y and earnings bin z) and option bw not be specified (pre-binned data)."
							exit 301
						}

						/* savebins() -- write the histogram + raw bin co-visitation
						   matrix to a dataset that vce(cluster stub) can read back.
						   Only meaningful when we are computing that matrix, i.e.
						   vce(cluster clustvar) on individual data. */
						local savebinsfile ""
						local savebinsrepl ""
						if `"`savebins'"' != "" {
							if `nvars'!=1 | "`clustvar'"=="" {
								noi di as error "savebins() requires vce(cluster clustvar) with individual-level data."
								exit 198
							}
							gettoken savebinsfile _sbrest : savebins, parse(",")
							local savebinsfile = strtrim(`"`savebinsfile'"')
							if `"`_sbrest'"' != "" {
								gettoken _sbc _sbrest : _sbrest, parse(",")
								local _sbrest = strtrim(`"`_sbrest'"')
								if `"`_sbrest'"' != "" {
									if `"`_sbrest'"' != "replace" {
										noi di as error "savebins() suboption must be replace."
										exit 198
									}
									local savebinsrepl replace
								}
							}
							if `"`savebinsfile'"' == "" {
								noi di as error "savebins() requires a file name."
								exit 198
							}
						}
						if `nvars'==2 { //find bw in pre-binned data
							loc y: word 1 of `varlist'
							loc z: word 2 of `varlist'
							
							sort `z'
							tempvar tmp
							gen `tmp'=`z'-`z'[_n-1]
							su `tmp'
							if r(Var)>0.01*r(mean) {
								noi di as error "Bandwidth differs in pre-binned data"
								exit 301
							}
							else loc bw=r(mean)
							sum `y'
							loc N=r(sum)
						}
						else { // collapse data
							loc z `varlist'
							tempvar bin y binid u

							if "`drop'"!="nodrop" {
								su `z'
								loc min=r(min)
								loc max=r(max)
								if abs(`cutoff'-floor((`z'-`min')/`bw')*`bw'-`min')>`bw'/10 {
									drop if `z'<`cutoff'-floor((`cutoff'-`min')/`bw')*`bw'
								}
								if abs(`cutoff'+floor((`max'-`cutoff')/`bw')*`bw'-`max')>`bw'/10 {
									drop if `z'>`cutoff'+floor((`max'-`cutoff')/`bw')*`bw'
								}
							}

							count
							
							loc N = r(N)
							tempvar u
							gen double `u' = (`z' - `cutoff') / `bw'

							/*
								Convention:
								  z == cutoff belongs to binid 0,
								  with center cutoff - bw/2.
								This makes the default limits(1 0) exclude the bunching bin.
							*/
							gen long `binid' = ceil(`u')
							replace `binid' = 0 if abs(`z' - `cutoff') < max(1e-12, abs(`bw')*1e-10)

							/* stash the individual (cluster, bin) pairs before the
							   collapse destroys them -- used to build the bin
							   co-visitation matrix for vce(cluster) */
							if "`clustvar'"!="" {
								mata: pbx_cluster_stash("`clustvar'", "`binid'")
							}

							collapse (count) `y'=`z', by(`binid')

							if "`zero'" != "nozero" {
								quietly summarize `binid', meanonly
								local idmin = r(min)
								local idmax = r(max)

								tempfile collapsed
								save `collapsed', replace

								clear
								set obs `=`idmax' - `idmin' + 1'
								gen long `binid' = `idmin' + _n - 1

								merge 1:1 `binid' using `collapsed', nogen
								replace `y' = 0 if missing(`y')
								sort `binid'
							}

							gen double `z' = `cutoff' + (`binid' - 0.5)*`bw'
						}
					

						if "`limits'" == "" {
							local L = 1
							local H = 0
						}
						else {
							gettoken L H : limits
							local L = real("`L'")
							local H = real("`H'")
						}

						tempvar zleft zright relbin crossbin edgehit inbunch

						gen double `zleft'  = `z' - `bw'/2
						gen double `zright' = `z' + `bw'/2

						gen byte `edgehit' = ///
							abs(`zleft'  - `cutoff') < 1e-8 | ///
							abs(`zright' - `cutoff') < 1e-8

						quietly count if `edgehit'
						local cutoff_on_edge = (r(N) > 0)

						gen int `relbin' = .
						gen byte `crossbin' = (`zleft' < `cutoff' & `zright' > `cutoff')
						
						if `cutoff_on_edge' {
							/*
								cutoff is a bin edge:
									relbin = -1 closest bin below cutoff
									relbin =  1 closest bin above cutoff
							*/

							replace `relbin' = -round((`cutoff' - `zright')/`bw') - 1 ///
								if `zright' <= `cutoff' + 1e-8

							replace `relbin' =  round((`zleft' - `cutoff')/`bw') + 1 ///
								if `zleft' >= `cutoff' - 1e-8

							gen byte `inbunch' = ///
								inrange(`relbin', -`L', -1) | ///
								inrange(`relbin',  1,  `H')
						}
						else {
							/*
								cutoff lies inside exactly one bin:
									relbin = 0 cutoff-containing bin
									relbin = -1 closest whole bin below it
									relbin =  1 closest whole bin above it

								Exclude cutoff bin plus L whole bins below and H whole bins above.
							*/
							quietly count if `crossbin'
							if r(N) != 1 {
								noi di as error "Expected exactly one bin crossing cutoff(). Check bw(), cutoff(), and bin construction."
								exit 301
							}

							quietly summarize `z' if `crossbin', meanonly
							local zcross = r(mean)

							replace `relbin' = round((`z' - `zcross')/`bw')

							gen byte `inbunch' = inrange(`relbin', -`L', `H')
						}

						/*
							Actual excluded-region edges from selected whole bins.
						*/
						quietly summarize `zleft' if `inbunch', meanonly
						if r(N) == 0 {
							noi di as error "No bins in the excluded region. Check limits(), cutoff(), and bw()."
							exit 301
						}
						local zL_excl_orig = r(min)

						quietly summarize `zright' if `inbunch', meanonly
						local zH_excl_orig = r(max)

						/*
							Bunch indicator.
						*/
						tempvar bunch
						egen `bunch' = group(`z') if `inbunch'
						replace `bunch' = 0 if missing(`bunch')
						quietly levelsof `bunch' if `bunch' > 0, local(bunchlevels)
						local Nbunch : word count `bunchlevels'

						//NORMALIZE Z
						tempvar z_orig
						gen double `z_orig'=`z' 
						loc cutoff_orig = `cutoff'
						loc bw_orig = `bw'
						
						if "`normalize'" != "nonormalize" {
							su `z'
							loc zmin_est=r(min)
							loc zmax_est=r(max)
							local zmid   = (`zmin_est' + `zmax_est')/2
							local xscale = (`zmax_est' - `zmin_est')/2

							replace `z' = (`z' - `zmid') / `xscale'

							local cutoff_est = (`cutoff' - `zmid') / `xscale'
							local bw_est     = `bw' / `xscale'
							
							local zL_excl_est = (`zL_excl_orig' - `zmid') / `xscale'
							local zH_excl_est = (`zH_excl_orig' - `zmid') / `xscale'
						}
						else {
							local cutoff_est = `cutoff'
							local bw_est = `bw'
							local xscale = 1
							local zL_excl_est = `zL_excl_orig'
							local zH_excl_est = `zH_excl_orig'
						}
						
						tempname table
						mkmat `y' `z', matrix(`table')
						mat colnames `table'= freq `z'

						/* raw (un-normalized) histogram: bin count + bin midpoint
						   in the original z units, exactly as supplied.  Used by
						   -polbunch_contrast- to refit stored specifications on
						   resampled histograms without a lossy normalize round
						   trip. */
						tempname binsraw
						mkmat `y' `z_orig', matrix(`binsraw')
						mat colnames `binsraw' = freq midpoint

						summarize `z', meanonly
						local zbar_est = r(max) + 0.5*`bw_est'
						
						tempvar side
						gen byte `side' = .

						local tol = max(1e-8, abs(`bw_orig')*1e-8)

						replace `side' = -1 if `bunch' == 0 & `zright' <= `zL_excl_orig' + `tol'
						replace `side' =  1 if `bunch' == 0 & `zleft'  >= `zH_excl_orig' - `tol'

						count if `bunch' == 0 & missing(`side')
						if r(N) > 0 {
							noi di as error "Some non-excluded bins cannot be classified as left or right of excluded region."
							noi list `z' `zleft' `zright' `relbin' `bunch' if `bunch' == 0 & missing(`side'), noobs
							exit 301
						}
						count if `side' == -1
						if r(N) == 0 {
							noi di as error "No bins below the excluded region."
							exit 301
						}

						count if `side' == 1
						if r(N) == 0 & `estimator' > 0 {
							noi di as error "No bins above the excluded region."
							exit 301
						}

						/*
							vce(cluster): build M, the bin co-visitation matrix.

							  _brow : each bin's position 1..J in bin order
							  _srow : each bin's row in the stacked estimating system
							          -- reference bins keep their own row (in bin
							          order, matching make_ystack), excluded bins all
							          fold into the trailing bunching-mass row.

							The raw J x J matrix M_raw (over all bins) comes either
							from the stashed individual visits (individual data) or
							from the user-supplied co-visitation stub (pre-binned
							data); pbx_covis_collapse then folds it onto the stacked
							system.  Everything downstream is identical for the two.
						*/
						if "`clustvar'"!="" {
							tempvar _isref _srow _brow
							sort `z'
							gen long `_brow'  = _n
							gen byte `_isref' = (`bunch'==0)
							gen long `_srow'  = sum(`_isref')
							count if `_isref'
							local _nref   = r(N)
							local _jstack = `_nref' + 1
							replace `_srow' = `_jstack' if !`_isref'
							local _nbin = _N
							tempname _srowmat _Mraw _Mstack

							if `nvars'==1 {
								tempname _bpos
								mkmat `binid' `_brow', matrix(`_bpos')
								mata: pbx_covis_build("`_bpos'", `_nbin')
								matrix `_Mraw' = r_pbx_Mraw
								local _nclust = r_pbx_nclust
							}
							else {
								/* resolve the co-visitation stub -> J variables
								   `clustvar'1..`clustvar'`_nbin', row-aligned with
								   the histogram (a full J x J matrix; the
								   diagonal-only shortcut is not offered because
								   diag(M) - y y'/n is not a valid covariance --
								   use vce(robust) as the histogram-only hedge). */
								local _mvars ""
								forvalues _k = 1/`_nbin' {
									capture confirm variable `clustvar'`_k'
									if _rc {
										noi di as error "vce(cluster `clustvar') on pre-binned data needs `_nbin' variables `clustvar'1..`clustvar'`_nbin' holding the J x J bin co-visitation matrix, row-aligned with the histogram. (`clustvar'`_k' not found.)"
										noi di as error "If you only have per-bin repeat counts, use vce(robust) on the histogram instead."
										exit 198
									}
									local _mvars `_mvars' `clustvar'`_k'
								}
								mkmat `_mvars', matrix(`_Mraw')
								mata: pbx_covis_check("`_Mraw'", "`y'")
								if r_pbx_covchk == 1 {
									noi di as error "vce(cluster `clustvar'): the supplied co-visitation matrix is not symmetric."
									exit 198
								}
								if r_pbx_covchk == 2 {
									noi di as error "vce(cluster `clustvar'): diag(M) has an entry below the corresponding bin count -- M[j,j] = sum_i c_ij^2 >= y_j must hold."
									exit 198
								}
								if r_pbx_covchk == 3 {
									noi di as error "vce(cluster `clustvar'): the co-visitation matrix is not `_nbin' x `_nbin'."
									exit 198
								}
								local _nclust = `nclusters'
								if `_nclust' > `N' {
									noi di as error "nclusters() (`_nclust') exceeds the total count (`N')."
									exit 198
								}
							}

							if "`savebins'"!="" {
								tempname _sbhist
								mkmat `y' `z_orig', matrix(`_sbhist')
								matrix colnames `_sbhist' = freq midpoint
							}

							mkmat `_srow', matrix(`_srowmat')
							mata: pbx_covis_collapse("`_Mraw'", "`_srowmat'", `_jstack')
							matrix `_Mstack' = r_pbx_Mstack

							if `_nclust' < 2 {
								noi di as error "vce(cluster): fewer than 2 clusters."
								exit 198
							}
						}
						local _covopt ""
						if "`clustvar'"!="" local _covopt covis(`_Mstack') nclust(`_nclust')

						//gen dummies
						tempvar dum dum2
						gen byte `dum' = `z' > `cutoff_est'
						gen byte `dum2' = `dum'
						count
						loc numbins=r(N)
						

						//Evaluate multicollinearity & estimate unrestricted model
							local rhsvars
							local coleq0
							local coleq1

							if `polynomial' > 0 {
								forvalues i = 1/`polynomial' {
									if `i' == 1 local rhsvars c.`z'
									else local rhsvars `rhsvars'##c.`z'

									local coleq0 `coleq0' h0
									local coleq1 `coleq1' h1
								}

								fvexpand `rhsvars'
								local names `r(varlist)' _cons
							}
							else {
								local rhsvars
								local names _cons
							}

							local coleq0 `coleq0' h0
							local coleq1 `coleq1' h1


							if `polynomial' > 0 {
								local nmiss = 1

								while `nmiss' {
									regress `y' 0.`dum'#(`rhsvars') 0.`dum' 1.`dum2'#(`rhsvars') 1.`dum2' if `bunch' == 0, nocons

									/*
										norankcheck: skip the separate-sides
										identification check entirely and keep the
										requested polynomial order.  A restricted
										estimator (2/3) can be identified at an order
										where two free one-sided polynomials are not;
										the unrestricted fit above is still run once
										to seed the profile starting values.
										Specification tests are disabled further
										below because they DO need an identified
										unrestricted fit.  The analytical bias does
										NOT -- polbunchbias reads the counterfactual
										polynomial, tax rates, window and elasticity
										from e() and never touches the separate
										one-sided polynomials -- so it is still run.
									*/
									if "`rankcheck'" == "norankcheck" {
										local nmiss = 0
										local rankforced 1
									}
									else {
									local nmiss = e(rank) < (`polynomial' + 1)*2

									if `nmiss' {
										if "`rankred'" != "norankred" {
											local note note
											local polynomial = `polynomial' - 1

											if `polynomial' < 0 {
												noi di as err "Could not estimate separate polynomials on either side of the cutoff."
												exit 301
											}

											// rebuild RHS and names
											local rhsvars
											local coleq0
											local coleq1

											if `polynomial' > 0 {
												forvalues i = 1/`polynomial' {
													if `i' == 1 local rhsvars c.`z'
													else local rhsvars `rhsvars'##c.`z'

													local coleq0 `coleq0' h0
													local coleq1 `coleq1' h1
												}

												fvexpand `rhsvars'
												local names `r(varlist)' _cons
											}
											else {
												local rhsvars
												local names _cons
											}

											local coleq0 `coleq0' h0
											local coleq1 `coleq1' h1
										}
										else {
											/*
												norankred: the two free one-sided
												polynomials are rank-deficient at this
												order, but the user asked to keep it.
												The restricted profile (est 2/3) can
												still be identified, but the
												minimum-distance / Hausman tests
												compare against this unrestricted fit
												-- which is now degenerate -- so they
												are disabled, exactly as under
												norankcheck.
											*/
											local nmiss = 0
											local rankforced_red 1
										}
									}
									} // end else (rank check performed)
								}
							}
								if "`note'"=="note" {
									noi di as text "Note: Polynomial order lowered to `polynomial' because of multicollinearity problems with the specified polynomial."
								}

								if "`rankforced'" == "1" {
									noi di as text "Note: norankcheck -- estimating with the requested polynomial(`polynomial') and no"
									noi di as text "      separate-sides rank check. Specification tests are disabled (they require"
									noi di as text "      an identified unrestricted fit); the analytical bias is still reported."
									local test none
								}

								if "`rankforced_red'" == "1" {
									noi di as text "Note: norankred -- the two one-sided polynomials are rank-deficient at polynomial(`polynomial'),"
									noi di as text "      but the order was kept as requested. Specification tests are disabled (they compare"
									noi di as text "      against the now-degenerate unrestricted fit); the restricted estimate and analytical"
									noi di as text "      bias are still reported."
									local test none
								}

						if inlist(`estimator',2,3) { //STARTING VALUES
							tempname h0coefs h1coefs bu
							mat `bu'=e(b)
							local Kb = `polynomial' + 1

							matrix `h0coefs' = `bu'[1, 1..`Kb']
							matrix `h1coefs' = `bu'[1, `=`Kb'+1'..`=2*`Kb'']
										/*
							Simple delta initializer from one beta/gamma coefficient relation.
							h0coef/h1coef ordering:
								beta_1 ... beta_K beta_0
						*/

						scalar dstart = .
						
						if "`positive'" != "" {
							local dstart_lower = 0
							local dstart_default = 0.05
							}
						else {
							local dstart_lower = -0.99
							local dstart_default = 0
							}

						if `estimator' == 2 {
							/*
								Prefer constant if usable, otherwise first polynomial coefficient.
							*/
							local j = `=`polynomial' + 1'

							if abs(`h1coefs'[1,`j']) > 1e-12 {
								scalar dstart = `h0coefs'[1,`j'] / `h1coefs'[1,`j'] - 1
							}

							if missing(dstart) | dstart <= 0 {
								local j = 1
								if `polynomial' >= 1 & abs(`h1coefs'[1,`j']) > 1e-12 {
									scalar dstart = `h0coefs'[1,`j'] / `h1coefs'[1,`j'] - 1
								}
							}
						}
						else if `estimator' == 3 & "`log'" == "" {
							/*
								Level case:
									gamma_j = beta_j * (1+delta)^(j+1)

								Use first usable polynomial coefficient, not the constant.
							*/
							forvalues j = 1/`polynomial' {
								if missing(dstart) | dstart <= 0 {
									if abs(`h0coefs'[1,`j']) > 1e-12 & ///
									   `h1coefs'[1,`j'] / `h0coefs'[1,`j'] > 0 {
										scalar dstart = ///
											(`h1coefs'[1,`j'] / `h0coefs'[1,`j'])^(1/(`j' + 1)) - 1
									}
								}
							}
						}
						else if `estimator' == 3 & "`log'" == "log" {
							/*
								Log case:
									gamma_0 - beta_0 ~= beta_1 * ln(1+delta)

								Requires polynomial >= 1 and beta_1 nonzero.
							*/
							if `polynomial' >= 1 & abs(`h0coefs'[1,1]) > 1e-12 {
								scalar dstart = exp( ///
									(`h1coefs'[1,`=`polynomial' + 1'] - ///
									 `h0coefs'[1,`=`polynomial' + 1']) / ///
									 `h0coefs'[1,1] ///
								) - 1
							}
						}

						/*
							Fallback / bounds.
						*/
						if missing(dstart) | dstart <= `dstart_lower' {
							scalar dstart = 0.05
						}

						if dstart > 0.5 {
							scalar dstart = 0.5
						}
							
						local dstart = scalar(dstart)
						} 
						else local dstart=0
					
							
						//BOOTSTRAP SETUP
						loc _wildflag = ("`boottype'"=="wild")
						loc _wildcode = cond("`wildwt'"=="webb",3,cond("`wildwt'"=="mammen",2,1))
						if "`vce'"=="bootstrap" {
							tempname p yorig
							if "`boottype'"=="bayesian" gen double `p'=`y'/`N'
							else if "`boottype'"=="multinomial" {
								gen double `yorig'=`y'
								recast double `y'
							}
							else {
								/* residual / wild bootstrap over bins.  The
								   reference distribution -- the fitted mean
								   m(theta-hat) that residuals are drawn around
								   and added back to -- is the SAME model as the
								   point estimate (h0 below the kink, h1 above),
								   not a separate per-side OLS.  It is not
								   available until theta-hat is in hand, so the
								   pool is built right after the s==0 fit below;
								   here we only widen `y' to double. */
								recast double `y'
								capture matrix drop r_mu_bins
								capture scalar drop r_rbsetup_ok
							}
						}
								
						local dotest = inlist(`estimator', 1, 2, 3,4) & "`test'" != "none" & "`vce'"!="none"
						
						//ESTIMATION AND INFERENCE
						tempname b V bs bmain Vmain b0 V0 b0s GR_raw y_raw GU_raw muU_raw GbcR_raw GbcU_raw
					tempname _crefs _bootflags

						/* failed-replication accounting (bootstrap vce); populated in the loop / summary below */
						local _nfail_rank    = 0
						local _nfail_draw    = 0
						local _nboot_ok      = .
						local _nfail_unrestr = .
						
						if inlist("`vce'","none","analytic") loc stop=0
						else loc stop=`bootreps'
						forvalues s=0/`stop' {
							if `s'==1 & "`nodots'"=="" nois _dots 0, title("Performing bootstrap repetitions...") reps(`bootreps')
							if `s'>0 { //resample outcome
								if "`boottype'"=="bayesian" {
									loc i=0
									loc factor=0
									loc obs=`N'
									while `obs'>0&`i'<`=_N-1' {
										loc ++i
										replace `y'=rbinomial(`obs',`p'/(1-`factor')) in `i'
										loc factor=`factor'+`p'[`i']
										loc obs=`obs'-`y'[`i']
									}
									if `obs'>0 replace `y'=`obs' in `=_N'
									else replace `y'=0 if _n>`i'
									if "`zero'"=="nozero" replace `y'=. if `y'==0
									}
								else if "`boottype'"=="multinomial" {
									replace `y'=rgamma(`yorig',1)
									su `y'
									replace `y'=`y'*`N'/r(sum)
								}
								else {
									mata: pbx_resboot_draw("`y'", `_wildflag', `_wildcode')
								}

							}
							
							//estimate model: single stacked profile branch for estimators 0/1/2/3

								/*
									vce(bootstrap, residual|wild) needs the fitted
									mean m(theta-hat) and its dispersion from the
									MAIN fit to build the residual pool.  Ask the
									s==0 fit for the analytic machinery (it fills
									r_mu_bins / e(dispersion); the analytic V it
									also computes is discarded -- `vce' stays
									"bootstrap").  Replications (s>0) never need it.
								*/
								local _vce_call vce(`vce')
								if `s'==0 & "`vce'"=="bootstrap" & inlist("`boottype'","residual","wild") ///
									local _vce_call vce(analytic)

								/*
									FAILED-REPLICATION ACCOUNTING (bootstrap vce).
									Separate-sides polynomial identification on THIS
									resample -- the same rank check polbunch runs once
									on the observed histogram (polbunch.ado, "Evaluate
									multicollinearity").  Inside the bootstrap it is
									diagnostic only: qrsolve() still returns a
									(pseudo-inverse) point estimate, so a rank-deficient
									draw need not surface as a missing coefficient -- but
									the two free one-sided polynomials, and hence the
									unrestricted counterfactual and every specification
									test built on it, are not identified on that draw.
									It never lowers the order or drops the draw here.
								*/
								local _rankfail = 0
								if `s' > 0 & `polynomial' > 0 & `estimator' != 4 {
									capture quietly regress `y' 0.`dum'#(`rhsvars') 0.`dum' 1.`dum2'#(`rhsvars') 1.`dum2' if `bunch' == 0, nocons
									if _rc | (e(rank) < (`polynomial' + 1)*2) local _rankfail = 1
								}
								if `s' > 0 & `bootreps' > 1 matrix `_bootflags' = nullmat(`_bootflags') \ (`_rankfail')

								if `estimator' == 4 {
									bunch_saez `y' `z' `side' `bunch', ///
										cutoff_orig(`cutoff_orig') ///
										bw_orig(`bw_orig') ///
										zl_excl_orig(`zL_excl_orig') ///
										zh_excl_orig(`zH_excl_orig') ///
										`_vce_call' hctype(`hctype') masscorr(`masscorr') `_covopt'
								}
								else {
									bunch_profile `y' `z' `side' `bunch', ///
										estimator(`estimator') k(`polynomial') ///
										cutoff_orig(`cutoff_orig') bw_orig(`bw_orig') ///
										cutoff_est(`cutoff_est') bw_est(`bw_est') ///
										l(`L') h(`H') ///
										`log' `normalize' `_vce_call' hctype(`hctype') masscorr(`masscorr') `_covopt' ///
										initdelta(`dstart') deltamax(`deltamax') ///
										zbar_est(`zbar_est') ///
										zl_excl_orig(`zL_excl_orig') ///
										zh_excl_orig(`zH_excl_orig') ///
										zl_excl_est(`zL_excl_est') ///
										zh_excl_est(`zH_excl_est') ///
										nosplit(`nosplit') ///
										`positive'
								}

								if `s'==0 & "`vce'"=="bootstrap" & inlist("`boottype'","residual","wild") {
									/* build the residual pool from m(theta-hat) */
									local _rbphi = .
									capture confirm scalar e(dispersion)
									if !_rc local _rbphi = e(dispersion)
									capture confirm matrix r_mu_bins
									if _rc {
										mata: pbx_resboot_setup("`y'","`z'","`side'","`bunch'", `polynomial')
									}
									else {
										capture scalar drop r_rbsetup_ok
										mata: pbx_resboot_setup_model("`y'", "`bunch'", `_rbphi')
										capture confirm scalar r_rbsetup_ok
										if _rc | r_rbsetup_ok!=1 {
											di as text "Note: vce(bootstrap, `boottype') -- the fitted counterfactual was" ///
												" unavailable for the residual pool; falling back to a per-side polynomial reference fit."
											mata: pbx_resboot_setup("`y'","`z'","`side'","`bunch'", `polynomial')
										}
									}
								}

								if `s'==0 {
									local phi_main = .
									capture confirm scalar e(dispersion)
									if !_rc local phi_main = e(dispersion)

									/* reference-region goodness of fit from the
									   MAIN fit -- captured before bunch_transform
									   / the outer ereturn re-post wipe it. */
									local _goflist gof_nbins gof_np gof_df deviance deviance_p ///
										pearson_x2 pearson_x2_p gof_ll deviance_null gof_ll_null ///
										r2_dev aic bic qaic qaicc gof_rmse gof_massresid ///
										dispersion_below dispersion_above deviance_below deviance_above ///
										gof_df_below gof_df_above deviance_below_p deviance_above_p
									foreach _g of local _goflist {
										local gof_`_g' = .
										capture confirm scalar e(`_g')
										if !_rc local gof_`_g' = e(`_g')
									}

									/*
										estimator 2/3 profile: was the structural
										delta weakly identified (boundary of the
										deltamax() search or a multi-modal profile
										SSR)?  Captured from the MAIN fit and reused
										for every bunch_transform call so the posted
										b has a constant width across bootstrap reps.
									*/
									local weakid_main = 0
									capture confirm scalar e(delta_weakid)
									if !_rc local weakid_main = e(delta_weakid)
									local nbasin_main = .
									capture confirm scalar e(delta_nbasin)
									if !_rc local nbasin_main = e(delta_nbasin)

									/* the robust/hc2/hc3 sandwich is not guaranteed
									   PSD (see variance_robust); when it needed
									   clipping, Stata's own -ereturn post- would
									   otherwise have silently zeroed e(V) with only
									   a terse warning -- surface it here instead. */
									local psdclip_main = .
									capture confirm scalar e(vce_psdclip)
									if !_rc local psdclip_main = e(vce_psdclip)
									if `psdclip_main'==1 {
										noi di as text "Note: the robust/HC sandwich was not positive semi-definite" ///
											" and was projected onto the nearest valid (PSD) covariance" ///
											" before reporting standard errors -- see e(vce_psdclip)."
									}
								}

							//STORE RESTRICTED STACKED GRADIENT for the scalar Hausman test.
							if `dotest' & `s'==0 & "`vce'"=="analytic" & ///
							   inlist("`test'","hausman","all","forceall") & inlist(`estimator',1,2,3) {
								matrix `GR_raw' = e(G_stack)
								matrix `y_raw'  = e(y_stack)
							}
							
							////TRANSFORM ESTIMATES
							if "`transform'"!="notransform" {
								summarize `y' if `bunch' > 0, meanonly
								local Hstar_obs = r(sum)
								local taxopts
								if "`t0'" != "" & "`t1'" != "" {
									local taxopts t0(`t0') t1(`t1')
								}
								if inlist(`estimator',0,1,2,3) {
									bunch_transform `z', ///
										estimator(`estimator') ///
										k(`polynomial') ///
										cutofforig(`cutoff_orig') ///
										cutoffest(`cutoff_est') ///
										bworig(`bw_orig') ///
										bwest(`bw_est') ///
										xscale(`xscale') ///
										low(`L') ///
										high(`H') ///
										zlexcl(`zL_excl_orig') ///
										zhexcl(`zH_excl_orig') ///
										zlexclest(`zL_excl_est') ///
										zhexclest(`zH_excl_est') ///
										`log' ///
										`constant' ///
										`taxopts' ///
										`grad' ///
										`normalize' ///
										zbar(`zbar_est') ///
										massobs(`Hstar_obs') ///
										nosplit(`nosplit') ///
										weakid(`weakid_main')
								}
								else {
									saez_transform, zstarorig(`cutoff_orig') bworig(`bw_orig') t0(`t0') t1(`t1') `log' `grad' `constant'
									}
								}
								
								matrix `b' = e(b)

								if `s' == 0 {
									matrix `bmain' = `b'
									if "`vce'"=="analytic" matrix `Vmain' = e(V)

									/*
										Scalar Hausman: the restricted model's
										elasticity and the delta-method Jacobian
										e(G) of the transform (its last row is
										d elasticity / d raw stacked coefs).
									*/
									local _shaus_eR = .
									if `dotest' & "`vce'"=="analytic" & ///
									   inlist("`test'","hausman","all","forceall") & inlist(`estimator',1,2,3) {
										capture matrix `GbcR_raw' = e(G)
										capture local _shaus_eR = _b[bunching:elasticity]
										if _rc local _shaus_eR = .
									}

									/*
										bunch_transform sets e(hasresp), but
										it's about to be overwritten by
										polbunch's own final ereturn block
										below (and, for bootstrap/bayes vce,
										by later replications' calls to this
										same transform step) -- capture it
										from the MAIN (s==0) estimate only, so
										it can be forwarded once the loop
										finishes. saez_transform doesn't set
										e(hasresp) (it signals an unsolved
										draw differently, via an all-missing
										b/V rather than a narrower one), so
										default to missing rather than
										erroring when it's absent.
									*/
									capture confirm scalar e(hasresp)
									if !_rc local hasresp_main = e(hasresp)
									else    local hasresp_main = .
								}
								else if `bootreps'>1 {
									mat `bs'=nullmat(`bs') \ `b'
								}

								/*
									CONTRAST: refit the model-consistent efficient
									reference -- estimator 3, exact inversion,
									splitmass -- on THIS SAME resampled histogram
									(`y' is already the draw for replication `s'),
									and record its elasticity.  Paired with the
									reported estimator's elasticity draw (which is
									already accumulated in `bs' / `bmain'), the two
									columns give a bootstrap estimate of
									Var(e_chosen - e_ref) that carries their
									dependence.  A draw on which the reference
									profile is weakly identified or the exact
									inversion has no real root contributes a
									missing value and is dropped pairwise.
								*/
								if `docontrast' {
									local _cref = .
									capture bunch_profile `y' `z' `side' `bunch', ///
										estimator(3) k(`polynomial') ///
										cutoff_orig(`cutoff_orig') bw_orig(`bw_orig') ///
										cutoff_est(`cutoff_est') bw_est(`bw_est') ///
										l(`L') h(`H') ///
										`log' `normalize' vce(none) hctype(-1) masscorr(1) ///
										initdelta(`dstart') deltamax(`deltamax') ///
										zbar_est(`zbar_est') ///
										zl_excl_orig(`zL_excl_orig') ///
										zh_excl_orig(`zH_excl_orig') ///
										zl_excl_est(`zL_excl_est') ///
										zh_excl_est(`zH_excl_est') ///
										nosplit(0) `positive'
									if !_rc {
										local _cwk = 0
										capture confirm scalar e(delta_weakid)
										if !_rc local _cwk = e(delta_weakid)
										if !`_cwk' {
											summarize `y' if `bunch' > 0, meanonly
											local _cHobs = r(sum)
											capture bunch_transform `z', ///
												estimator(3) k(`polynomial') ///
												cutofforig(`cutoff_orig') cutoffest(`cutoff_est') ///
												bworig(`bw_orig') bwest(`bw_est') xscale(`xscale') ///
												low(`L') high(`H') ///
												zlexcl(`zL_excl_orig') zhexcl(`zH_excl_orig') ///
												zlexclest(`zL_excl_est') zhexclest(`zH_excl_est') ///
												`log' `taxopts' `normalize' ///
												zbar(`zbar_est') massobs(`_cHobs') nosplit(0) weakid(0)
											if !_rc {
												capture local _cref = _b[bunching:elasticity]
												if _rc local _cref = .
											}
										}
									}
									if `s'==0 local _cref0 = `_cref'
									else if `bootreps'>1 mat `_crefs' = nullmat(`_crefs') \ (`_cref')
								}


								//IF TESTING: ALSO ESTIMATE UNRESTRICTED MODEL
								if `dotest'&`estimator'!=4 {
									bunch_profile `y' `z' `side' `bunch', ///
										estimator(0) k(`polynomial') ///
										cutoff_orig(`cutoff_orig') bw_orig(`bw_orig') ///
										cutoff_est(`cutoff_est') bw_est(`bw_est') ///
										l(`L') h(`H') ///
										`log' `normalize' vce(`vce') hctype(`hctype') masscorr(`masscorr') `_covopt' ///
										initdelta(`dstart') ///
										zbar_est(`zbar_est') ///
										zl_excl_orig(`zL_excl_orig') ///
										zh_excl_orig(`zH_excl_orig') ///
										zl_excl_est(`zL_excl_est') ///
										zh_excl_est(`zH_excl_est') ///
										`positive'

									if `s' == 0 {
										matrix `b0' = e(b)
										if "`vce'"=="analytic" matrix `V0' = e(V)
									}
									else if `bootreps' > 1 {
										matrix `b0s' = nullmat(`b0s') \ e(b)
									}

									/*
										Scalar Hausman: capture the unrestricted
										stacked gradient, then transform the
										estimator-0 fit to get its elasticity and
										the transform Jacobian.  Analytic vce only
										(the influence functions need e(G_stack));
										main fit only.  e() is left holding the
										transformed estimator-0 model, but nothing
										downstream reads it -- b0/V0 are already
										stored raw above.
									*/
									local _shaus_eU = .
									if `s'==0 & "`vce'"=="analytic" & ///
									   inlist("`test'","hausman","all","forceall") & inlist(`estimator',1,2,3) {
										capture matrix `GU_raw'  = e(G_stack)
										capture matrix `muU_raw' = e(mu_stack)
										summarize `y' if `bunch' > 0, meanonly
										local _Hobs_u = r(sum)

										/*
											e_U -- the elasticity the unrestricted
											two-sided fit implies -- always on the
											EXACT + SPLITMASS axes, the consistent
											combination (estimators 0 and 3 are
											unbiased there), whatever inversion/mass
											axes the reported estimate used.  So for
											the model-consistent estimator the contrast
											isolates the cross-kink density restriction;
											for the naive / Chetty estimators it also
											reflects their own constant/poolmass
											approximations, which is appropriate --
											those are part of the estimator's bias.  If
											the exact inversion has no real root on the
											extrapolated unrestricted counterfactual,
											e_U is left missing and the Hausman test is
											not reported (lower polynomial(), narrow
											limits(), or read the minimum-distance
											test, which needs no inversion).
										*/
										capture bunch_transform `z', ///
											estimator(0) k(`polynomial') ///
											cutofforig(`cutoff_orig') cutoffest(`cutoff_est') ///
											bworig(`bw_orig') bwest(`bw_est') xscale(`xscale') ///
											low(`L') high(`H') ///
											zlexcl(`zL_excl_orig') zhexcl(`zH_excl_orig') ///
											zlexclest(`zL_excl_est') zhexclest(`zH_excl_est') ///
											`log' `taxopts' `normalize' ///
											zbar(`zbar_est') massobs(`_Hobs_u') nosplit(0) weakid(0)
										if !_rc {
											capture matrix `GbcU_raw' = e(G)
											capture local _shaus_eU = _b[bunching:elasticity]
											if _rc local _shaus_eU = .
										}
									}
								}
								if `s' > 0 & "`nodots'"=="" noi _dots `s' 0
								
							}
							
								
							//bootstrap inference sunmmary & test
							if "`vce'"=="bootstrap" {
								/*
									clear+svmat+corr needs an empty dataset to
									load the bootstrap draws into, but the
									TEST RESTRICTIONS and POST RESULTS code
									below still needs the original estimation
									sample (z, side, ...) back afterward.
									-preserve- can't be used here: Stata only
									allows one preserve active at a time, and
									this whole command is already running
									inside its own outer preserve, so a second
									one errors with "already preserved". Round
									-trip through a tempfile instead, which
									doesn't touch the preserve stack -- without
									it, z stays cleared for the rest of the
									command, surfacing later as "variable z
									not found".
								*/
								tempfile _pbsim_boot_data
								quietly save `_pbsim_boot_data'

								clear
								svmat double `bs'
								corr _all, cov
								mat `Vmain'=r(C)

								/*
									FAILED-REPLICATION COUNT.
									  _nfail_draw -- replications whose REPORTED coefficient
									     vector came back missing on an element the point
									     estimate identifies (no real root in the response
									     inversion, a boundary / multi-modal delta, a hard
									     fit failure).  `corr _all, cov` above and the
									     percentile / BC quantiles below already drop these
									     pairwise; here we only count them.  e(bootreps_ok)
									     / e(bootreps_fail) are defined off this count.
									  _nfail_rank -- resamples on which the two free
									     one-sided polynomials are not separately identified
									     (from `_bootflags', set in the loop).  qrsolve()
									     still returns a minimum-norm fit, so for the naive
									     estimators these enter the (co)variance as
									     degenerate draws; for the restricted estimators
									     (2/3) only the specification test's unrestricted
									     fit is affected.  Reported separately, advisory.
								*/
								local _nk = colsof(`bmain')
								tempvar _bootbad _bootrk
								quietly gen byte `_bootbad' = 0
								forvalues j = 1/`_nk' {
									if !missing(`bmain'[1,`j']) quietly replace `_bootbad' = 1 if missing(`bs'`j')
								}
								quietly gen byte `_bootrk' = 0
								capture confirm matrix `_bootflags'
								if !_rc {
									local _nrf = rowsof(`_bootflags')
									forvalues _fs = 1/`_nrf' {
										if `_bootflags'[`_fs',1]==1 quietly replace `_bootrk' = 1 in `_fs'
									}
								}
								quietly count if `_bootbad'
								local _nfail_draw = r(N)
								quietly count if `_bootrk'
								local _nfail_rank = r(N)
								local _nboot_ok   = `bootreps' - `_nfail_draw'
								quietly drop `_bootbad' `_bootrk'

								/*
									Non-normal bootstrap CIs (bootci = percentile | bc).
									Stored in e(ci); the displayed table stays normal.
								*/
								if "`bootci'"!="normal" {
									local _nk = colsof(`bmain')
									tempname _ciM
									matrix `_ciM' = J(`_nk',2,.)
									forvalues j = 1/`_nk' {
										capture confirm variable `bs'`j'
										if _rc continue
										local _plo 2.5
										local _phe 97.5
										if "`bootci'"=="bc" {
											quietly count if !missing(`bs'`j')
											local _ntot = r(N)
											quietly count if `bs'`j' < `bmain'[1,`j'] & !missing(`bs'`j')
											local _frac = cond(`_ntot'>0, r(N)/`_ntot', 0.5)
											if `_frac'<=0 local _frac = 0.5/max(`_ntot',1)
											if `_frac'>=1 local _frac = 1 - 0.5/max(`_ntot',1)
											local _z0 = invnormal(`_frac')
											local _plo = 100*normal(2*`_z0' + invnormal(0.025))
											local _phe = 100*normal(2*`_z0' + invnormal(0.975))
										}
										quietly _pctile `bs'`j', p(`_plo' `_phe')
										matrix `_ciM'[`j',1] = r(r1)
										matrix `_ciM'[`j',2] = r(r2)
									}
									matrix colnames `_ciM' = ll ul
								}

								/*
									CONTRAST: combine the reported estimator's
									elasticity draws (column `_ecol' of `bs', matched
									to `bmain') with the paired reference draws in
									`_crefs' to form d* = e_chosen* - e_ref*.  The
									bootstrap SD of d* is the SE of the difference
									with the two estimators' dependence built in; the
									percentile interval and the bootstrap p-value
									(share of |d* - d0| >= |d0|) come from the same
									draws.  Missing on either side is dropped
									pairwise.
								*/
								if `docontrast' {
									local _c_ok = 0
									local _ecol = 0
									local _cj = 0
									if "`_cref0'"=="" local _cref0 = .
									local _cbmnm : colnames `bmain'
									foreach _nm of local _cbmnm {
										local ++_cj
										if "`_nm'"=="elasticity" local _ecol = `_cj'
									}
									capture confirm matrix `_crefs'
									local _chasref = (_rc==0)
									if `_ecol'>0 & !missing(`_cref0') & `_chasref' {
										capture confirm variable `bs'`_ecol'
										if !_rc {
											local _ce0 = `bmain'[1,`_ecol']
											local _d0  = `_ce0' - `_cref0'
											tempvar _cev _crv _cdv
											gen double `_cev' = `bs'`_ecol'
											gen double `_crv' = .
											local _crn = rowsof(`_crefs')
											forvalues _ci = 1/`_crn' {
												quietly replace `_crv' = `_crefs'[`_ci',1] in `_ci'
											}
											gen double `_cdv' = `_cev' - `_crv'
											quietly count if !missing(`_cdv')
											local _c_nused = r(N)
											if `_c_nused'>=2 & !missing(`_d0') {
												quietly summarize `_cdv'
												local _c_sed = r(sd)
												quietly correlate `_cev' `_crv' if !missing(`_cdv')
												local _c_corr = r(rho)
												quietly _pctile `_cdv', p(2.5 97.5)
												local _c_lo = r(r1)
												local _c_hi = r(r2)
												quietly count if !missing(`_cdv') & abs(`_cdv' - `_d0') >= abs(`_d0')
												local _c_ppct = r(N)/`_c_nused'
												if `_c_sed'>0 & !missing(`_c_sed') {
													local _c_z     = `_d0'/`_c_sed'
													local _c_pnorm = 2*normal(-abs(`_c_z'))
												}
												else {
													local _c_z     = .
													local _c_pnorm = .
												}
												local _c_ok = 1
											}
										}
									}
									if !`_c_ok' {
										if `_ecol'==0 ///
											local _c_why "the reported estimator's elasticity is not identified on this window (weak delta / no real root)"
										else if missing(`_cref0') ///
											local _c_why "the estimator(3) reference elasticity is not identified on the full sample"
										else ///
											local _c_why "fewer than 2 replications yielded a paired elasticity difference"
									}
								}

								quietly use `_pbsim_boot_data', clear

								if `dotest' {
									if `estimator'!=4 {
										clear
										svmat `b0s'
										corr _all, cov
										matrix `V0' = r(C)
										/* failed unrestricted (two-sided) draws feeding the bootstrap spec tests */
										local _nk0 = colsof(`b0')
										tempvar _bad0
										quietly gen byte `_bad0' = 0
										forvalues _fj = 1/`_nk0' {
											if !missing(`b0'[1,`_fj']) quietly replace `_bad0' = 1 if missing(`b0s'`_fj')
										}
										quietly count if `_bad0'
										local _nfail_unrestr = r(N)
										quietly use `_pbsim_boot_data', clear
									}

								}
							}
							

								
							//TEST RESTRICTIONS
			
							if `dotest' {
								if "`_shaus_eR'"=="" local _shaus_eR = .
								if "`_shaus_eU'"=="" local _shaus_eU = .
								local _shaus_seU = .

								if `estimator'==4 { //saez: Post main model to b0 V0
									mat `b0' = `bmain'
									mat `V0' = `Vmain'
								}

								/*
									Which tests to run.
									  wald            -- linear restriction (naive, Saez)
									  minimumdistance -- omnibus chi2_{K+1} on the cross-kink
									                     coefficient vector (Chetty, model-consistent)
									  hausman         -- focused chi2_1 on the elasticity
									For the naive estimator wald IS the minimum-distance
									test (no nuisance parameter), so it is not run twice.
								*/
								if inlist("`test'","all","forceall") {
									if `estimator'==4              local _tests_to_run "wald"
									else if `estimator'==1         local _tests_to_run "wald hausman"
									else if "`test'"=="forceall"   local _tests_to_run "wald minimumdistance hausman"
									else                           local _tests_to_run "minimumdistance hausman"
								}
								else {
									local _tests_to_run "`test'"
									local testname = cond("`test'"=="wald","Wald test", ///
										cond("`test'"=="minimumdistance","Minimum-distance test","Hausman test"))
								}

								/*
									The scalar Hausman test needs the analytic stacked
									gradients (e(G_stack) of both fits); drop it when
									vce() is not analytic.
								*/
								if "`vce'"!="analytic" & strpos(" `_tests_to_run' "," hausman ") {
									local _ttr2 ""
									foreach _tt of local _tests_to_run {
										if "`_tt'"!="hausman" local _ttr2 "`_ttr2' `_tt'"
									}
									local _tests_to_run = strtrim("`_ttr2'")
									if "`test'"=="hausman" {
										noi di as text "Note: the Hausman test requires vce(analytic|robust|hc2|hc3|cluster); not computed under vce(`vce')."
										local dotest = 0
									}
								}

								// Post b0/V0 once for all wald/md tests
								local _needs_post = 0
								foreach _tt of local _tests_to_run {
									if inlist("`_tt'","wald","minimumdistance") local _needs_post = 1
								}
								if `_needs_post' {
									local nm: colnames `b0'
									local neq: coleq `b0'
									mat colnames `V0'=`nm'
									mat rownames `V0'=`nm'
									mat coleq `V0'=`neq'
									mat roweq `V0'=`neq'
									ereturn post `b0' `V0'
									ereturn local properties "b V"
								}

								// Run each test and store per-test results
								local _tests_done ""
								foreach _tt of local _tests_to_run {

									if `estimator'!=4 {
										if "`_tt'"=="minimumdistance" {
											local init 0.05
											capture local init = _b[bunching:shift]
											if _rc | missing(real("`init'")) {
												capture local init = _b[bunching:delta]
											}
											if _rc | missing(real("`init'")) {
												local init 0.05
											}
											capture noisily polbunch_minimumdistancetest, ///
												estimator(`estimator') ///
												k(`polynomial') ///
												cutofforig(`cutoff_orig') ///
												cutoffest(`cutoff_est') ///
												bworig(`bw_orig') ///
												bwest(`bw_est') ///
												zbar(`zbar_est') ///
												`normalize' ///
												`log' ///
												`positive' ///
												initdelta(`init')
										}
										else if "`_tt'"=="wald" {
											capture noisily polbunch_waldtest, ///
												estimator(`estimator') ///
												k(`polynomial') ///
												cutofforig(`cutoff_orig') ///
												cutoffest(`cutoff_est') ///
												bworig(`bw_orig') ///
												bwest(`bw_est') ///
												zbar(`zbar_est') ///
												`normalize' ///
												`log'
										}
										else if "`_tt'"=="hausman" {
											capture noisily polbunch_shausmantest, ///
												gbcr(`GbcR_raw') gr(`GR_raw') er(`_shaus_eR') ///
												gbcu(`GbcU_raw') gu(`GU_raw') eu(`_shaus_eU') ///
												ystack(`y_raw') ///
												muu(`muU_raw') ///
												hctype(`hctype') ///
												`_covopt'
										}
									}
									else { //Saez: Simple Wald test of h0 vs h1
										if "`log'"=="" capture test _b[h1:_cons] - _b[h0:_cons] - (`bw_orig'/`cutoff_orig')*_b[bunching:number_bunchers] = 0
										else capture test _b[h1:_cons] = _b[h0:_cons]
									}

									local _test_rc = _rc
									local _fc = 0
									capture confirm scalar r(failcode)
									if !_rc {
										local _fc = r(failcode)
										if missing(`_fc') local _fc = 0
									}

									if (`_test_rc' != 0) | (`_fc' != 0) {
										local _tname = cond("`_tt'"=="wald","Wald test", ///
											cond("`_tt'"=="minimumdistance","Minimum-distance test","Hausman test"))
										if `_test_rc' != 0 local _tfc = `_test_rc'
										else local _tfc = `_fc'
										noi di as text "Note: `_tname' could not be computed; statistic not reported (rc=`_tfc')."
										if "`_tt'"=="hausman" & inlist(`_tfc',110,111) {
											noi di as text "      The exact response inversion has no real root on the unrestricted"
											noi di as text "      counterfactual.  Use a lower polynomial() or narrower limits(), or read"
											noi di as text "      the minimum-distance test, which needs no inversion."
										}
									}
									else {
										local chi2_`_tt' = r(chi2)
										local p_`_tt'    = r(p)
										local df_`_tt'   = r(df)
										if "`_tt'"=="minimumdistance" local delta_md = r(delta)
										if "`_tt'"=="hausman" local _shaus_seU = r(seU)
										local _tests_done "`_tests_done' `_tt'"
									}
								}

								local _tests_done = strtrim("`_tests_done'")
								if "`_tests_done'"=="" local dotest = 0
							}

						
						//POST RESULTS
						su `z' if `side'==-1, meanonly
						loc dL=r(mean)
						su `z' if `side'==1, meanonly
						loc dR=r(mean)

						su `z_orig', meanonly
						local _zlo_val = r(min) - `bw_orig'/2
						local _zhi_val = r(max) + `bw_orig'/2

						restore

						/*
							scale() -- glm-style overdispersion multiplier
							(vce(conventional) only; validated at parse time).
							scale(x2) uses the Pearson phi-hat from the reference
							bins; scale(#) uses the literal factor.
						*/
						local scaleapplied = 1
						if "`vce'"=="analytic" & `hctype'<0 & "`scalemode'"!="1" {
							if "`scalemode'"=="x2" {
								if !missing(`phi_main') & `phi_main'>0 local scaleapplied = `phi_main'
							}
							else local scaleapplied = real("`scalemode'")
							matrix `Vmain' = `scaleapplied' * `Vmain'
						}

						if "`vce'"!="none" {
							local nm: colnames `bmain'
							local neq: coleq `bmain'
							mat colnames `Vmain'=`nm'
							mat rownames `Vmain'=`nm'
							mat coleq `Vmain'=`neq'
							mat roweq `Vmain'=`neq'

							eret post `bmain' `Vmain', esample(`touse') depname(freq) obs(`N')
						}
						else eret post `bmain', esample(`touse') obs(`N') depname(freq)
						if `dotest' {
							foreach _tt of local _tests_done {
								estadd scalar chi2_`_tt' = `chi2_`_tt''
								estadd scalar p_`_tt'    = `p_`_tt''
								estadd scalar df_`_tt'   = `df_`_tt''
							}
							if strpos(" `_tests_done' "," minimumdistance ") estadd scalar delta_md = `delta_md'
							if strpos(" `_tests_done' "," hausman ") {
								capture confirm number `_shaus_eU'
								if !_rc & !missing(`_shaus_eU')  ereturn scalar elast_unrestricted    = `_shaus_eU'
								capture confirm number `_shaus_seU'
								if !_rc & !missing(`_shaus_seU') ereturn scalar se_elast_unrestricted = `_shaus_seU'
							}
						}
						if `docontrast' {
							ereturn local contrast_ref "estimator(3) exact splitmass"
							if `_c_ok' {
								ereturn scalar contrast_elast     = `_ce0'
								ereturn scalar contrast_elast_ref = `_cref0'
								ereturn scalar contrast_diff      = `_d0'
								ereturn scalar contrast_se        = `_c_sed'
								if !missing(`_c_z')     ereturn scalar contrast_z = `_c_z'
								if !missing(`_c_pnorm') ereturn scalar contrast_p = `_c_pnorm'
								ereturn scalar contrast_p_pctile  = `_c_ppct'
								ereturn scalar contrast_corr      = `_c_corr'
								ereturn scalar contrast_ci_ll     = `_c_lo'
								ereturn scalar contrast_ci_ul     = `_c_hi'
								ereturn scalar contrast_reps      = `bootreps'
								ereturn scalar contrast_reps_used = `_c_nused'
								ereturn scalar contrast_reps_fail = `bootreps' - `_c_nused'
							}
						}
						ereturn scalar polynomial=`polynomial'
						ereturn scalar lower_limit=`zL_excl_orig'
						ereturn scalar upper_limit=`zH_excl_orig'
						ereturn local normalize="`normalize'"
						ereturn scalar estimator=`estimator'
						/*
							Toggle axes, so a later stand-alone polbunchbias
							(e() mode) reproduces the same response-inversion
							and mass calculation this estimate used.
						*/
						ereturn scalar constant = `useconstant'
						ereturn scalar nosplit  = `nosplit'
						/*
							e(cmd) used to be set only when vce()!="none",
							which broke "estimates store"/"estimates
							restore" (and hence anything built on them,
							e.g. polbunchsim's vce(none)/btype(0) path)
							with "last estimation results not found" --
							there's no postestimation companion command
							here that would need e(V) to exist just
							because e(cmd) is set, so this is set
							unconditionally now.
						*/
						ereturn local cmd "polbunch"
						ereturn local cmdname "polbunch"
						ereturn local estat_cmd "polbunch_estat"
						ereturn local title 	"Polynomial bunching estimates"
						ereturn local cmdline 	`"polbunch `cmdline0'"'
						ereturn matrix table=`table'
						capture confirm matrix `binsraw'
						if !_rc ereturn matrix bins=`binsraw'
						ereturn local binname "`z'"
						ereturn scalar bw=`bw'
						ereturn scalar cutoff_orig = `cutoff_orig'
						ereturn scalar cutoff_est  = `cutoff_est'
						ereturn scalar bw_orig     = `bw_orig'
						ereturn scalar bw_est      = `bw_est'
						ereturn scalar xscale      = `xscale'
						ereturn scalar zL_excl_est = `zL_excl_est'
						ereturn scalar zH_excl_est = `zH_excl_est'
						ereturn scalar dL = `dL'
						ereturn scalar dR = `dR'
						ereturn local zname = "`z'"
						ereturn local transform="`transform'"

						/*
							e(vcetype) must stay short -- Stata prints it above
							"std. err." in the coefficient table.  The detail
							(wild/residual/hc flavour/scale) lives in e(vce),
							e(boottype), e(scale), etc.
						*/
						if "`clustvar'"!=""       estadd local vcetype "Robust"
						else if "`vce'"=="bootstrap"   estadd local vcetype "Bootstrap"
						else if `hctype'>=0       estadd local vcetype "Robust"
						else if "`scalemode'"!="1" estadd local vcetype "Scaled"
						else if "`transform'"=="notransform" estadd local vcetype "Analytic"
						else                      estadd local vcetype "Delta-method"
						/* inference-mode metadata */
						ereturn local vce "`vtype'"
						if "`clustvar'"!="" {
							ereturn scalar N_clust = `_nclust'
							if `nvars'==1 ereturn local clustvar "`clustvar'"
							else          ereturn local covisstub "`clustvar'"
						}
						if "`vce'"=="bootstrap" {
							ereturn local boottype "`boottype'"
							ereturn local bootci   "`bootci'"
							ereturn scalar bootreps = `bootreps'
							if "`boottype'"=="wild" ereturn local wildweights "`wildwt'"
							/* failed-replication accounting */
							capture confirm number `_nboot_ok'
							if !_rc & !missing(`_nboot_ok') {
								ereturn scalar bootreps_ok        = `_nboot_ok'
								ereturn scalar bootreps_fail      = `_nfail_draw'
								ereturn scalar bootreps_fail_rank = `_nfail_rank'
							}
							capture confirm number `_nfail_unrestr'
							if !_rc & !missing(`_nfail_unrestr') ereturn scalar bootreps_fail_unrestricted = `_nfail_unrestr'
							if "`bootci'"!="normal" {
								capture confirm matrix `_ciM'
								if !_rc {
									local _nm : colnames e(b)
									local _eq : coleq e(b)
									matrix rownames `_ciM' = `_nm'
									matrix roweq    `_ciM' = `_eq'
									ereturn matrix ci_`bootci' = `_ciM'
								}
							}
						}
						else if `hctype'>=0 {
							ereturn scalar masscorr = `masscorr'
							if !missing(`phi_main') ereturn scalar dispersion = `phi_main'
							capture confirm number `psdclip_main'
							if !_rc & !missing(`psdclip_main') ereturn scalar vce_psdclip = `psdclip_main'
						}
						else if "`vce'"=="analytic" {
							if !missing(`phi_main') ereturn scalar dispersion = `phi_main'
							if "`scalemode'"!="1" {
								ereturn local scale "`scalemode'"
								ereturn scalar scalefactor = `scaleapplied'
							}
						}

						/* reference-region goodness of fit, captured from the
						   s==0 fit.  Present whenever an analytic V was formed
						   -- vce(conventional|robust|hc*|cluster) and
						   vce(bootstrap, residual|wild) (whose s==0 draw runs
						   analytic); absent for vce(bootstrap, multinomial|
						   bayesian) and vce(none). */
						foreach _g of local _goflist {
							if !missing(`gof_`_g'') ereturn scalar `_g' = `gof_`_g''
						}

						if "`log'"=="log" ereturn scalar log=1
						else ereturn scalar log=0
						if "`t0'" != "" ereturn scalar t0 = `t0'
						if "`t1'" != "" ereturn scalar t1 = `t1'
						ereturn scalar zlo = `_zlo_val'
						ereturn scalar zhi = `_zhi_val'
						/*
							Forwards bunch_transform's e(hasresp) (captured
							from the main, s==0 estimate above, before it
							got overwritten) -- whether shift/marginal_response/
							elasticity are identified for this estimation.
							Missing when the transform step wasn't the
							bunch_transform path (e.g. Saez, estimator 4) or
							didn't run (notransform).
						*/
						capture confirm number `hasresp_main'
						if !_rc ereturn scalar hasresp = `hasresp_main'
						else     ereturn scalar hasresp = .

						/*
							estimator 2/3 only: delta_weakid = 1 when the
							structural-delta profile search hit the deltamax()
							boundary or found competing minima (multi-modal
							concentrated SSR).  Then hasresp is forced to 0 and
							shift/marginal_response/elasticity are withheld.
							delta_nbasin reports how many local minima the grid
							found (a diagnostic).
						*/
						if inlist(`estimator',2,3) {
							capture confirm number `weakid_main'
							if !_rc ereturn scalar delta_weakid = `weakid_main'
							capture confirm number `nbasin_main'
							if !_rc & !missing(`nbasin_main') ereturn scalar delta_nbasin = `nbasin_main'
							ereturn scalar delta_nonneg = ("`allownegative'" == "")
						}

						//Display results
						noi {
							di _newline
							di "`e(title)'"
							eret di
							if "`clustvar'"!="" {
								if `nvars'==1 di as txt "(Std. err. adjusted for `_nclust' clusters in `clustvar')"
								else          di as txt "(Std. err. adjusted for `_nclust' clusters; co-visitation matrix supplied)"
							}
							if "`vce'"=="bootstrap" & "`bootci'"!="normal" {
								local _cilab = cond("`bootci'"=="bc","bias-corrected","percentile")
								di as txt "Note: coefficient table shows normal-approximation CIs; `_cilab' bootstrap"
								di as txt "      CIs are stored in e(ci_`bootci')."
							}
							if "`vce'"=="bootstrap" {
								capture confirm number `_nfail_draw'
								if _rc local _nfail_draw = 0
								capture confirm number `_nfail_rank'
								if _rc local _nfail_rank = 0
								if (`_nfail_draw' > 0) | (`_nfail_rank' > 0) {
									di as txt "{hline 78}"
									if `_nfail_draw' > 0 {
										di as txt "Note: " as res "`_nfail_draw'" as txt " of `bootreps' replications returned a non-finite estimate (no real root in"
										di as txt "      the response inversion, a boundary/multi-modal delta, or a fit failure) and"
										di as txt "      were dropped pairwise from the (co)variance and the percentile/BC CIs."
										di as txt "      e(bootreps_ok) = " as res "`_nboot_ok'" as txt " of " as res "`bootreps'" as txt "."
									}
									if `_nfail_rank' > 0 {
										di as txt "Note: on " as res "`_nfail_rank'" as txt " of `bootreps' resamples the counterfactual polynomial could not be"
										di as txt "      estimated separately on both sides of the cutoff at polynomial(`polynomial')."
										if inlist(`estimator',0,1) {
											di as txt "      qrsolve() returns a minimum-norm fit for these, so they still enter the"
											di as txt "      bootstrap (co)variance -- consider a lower polynomial() or a pooled counterfactual."
										}
										else {
											di as txt "      The restricted estimate is unaffected; the specification test's unrestricted"
											di as txt "      two-sided fit is degenerate on those draws."
										}
										capture confirm number `_nfail_unrestr'
										if !_rc & `dotest' & !missing(`_nfail_unrestr') & `_nfail_unrestr' > 0 ///
											di as txt "      " as res "`_nfail_unrestr'" as txt " unrestricted-fit draws also came back non-finite and were dropped."
									}
									di as txt "{hline 78}"
								}
							}
							if !missing(`gof_r2_dev') {
								local _gr2 : di %5.3f `gof_r2_dev'
								local _gph : di %5.2f cond(`gof_gof_df'>0, `gof_pearson_x2'/`gof_gof_df', .)
								local _gqc : di %8.1f `gof_qaicc'
								di as txt "Counterfactual fit (reference bins):  deviance R2 " as res "`_gr2'" as txt "   phi-hat " as res "`_gph'" as txt "   QAICc " as res "`_gqc'"
								if !missing(`gof_dispersion_below') & !missing(`gof_dispersion_above') {
									local _gpb : di %5.2f `gof_dispersion_below'
									local _gpa : di %5.2f `gof_dispersion_above'
									di as txt "      phi-hat below the cutoff (A3 only) " as res "`_gpb'" as txt " | above (+ A1/A2) " as res "`_gpa'"
								}
								if inlist(`estimator',2,3) & !missing(`gof_gof_massresid') {
									if abs(`gof_gof_massresid') > 1e-4 {
										local _gmr : di %9.2e `gof_gof_massresid'
										di as txt "      (mass-restriction residual " as res "`_gmr'" as txt " relative -- delta solve not fully converged)"
									}
								}
								di as txt "      more: estat gof"
							}
							if `dotest' {
								tempname b
								matrix `b' = e(b)

								local stub = strlen("`e(depvar)'")

								local cn : colnames `b'
								local eq : coleq `b'

								foreach x of local cn {
									local stub = max(`stub', strlen("`x'"))
								}

								foreach x of local eq {
									local stub = max(`stub', strlen("`x'"))
								}

								local stub = max(`stub', 12)
								local W = `stub' + 67

								if inlist("`test'","all","forceall") {
									local _first = 1
									foreach _tt of local _tests_done {
										local _tshort = cond("`_tt'"=="wald","Wald", ///
											cond("`_tt'"=="minimumdistance","Min. dist.","Hausman"))
										if `_first' {
											di as txt "Model assumption tests:" ///
												_col(`=`W'-35') "`_tshort': Chi2(`df_`_tt'')" ///
												_col(`=`W'-10') as res %10.4f `chi2_`_tt''
											local _first = 0
										}
										else {
											di as txt _col(`=`W'-35') "`_tshort': Chi2(`df_`_tt'')" ///
												_col(`=`W'-10') as res %10.4f `chi2_`_tt''
										}
										di as txt _col(`=`W'-35') "p-value" ///
											_col(`=`W'-10') as res %10.4f `p_`_tt''
									}
								}
								else {
									di as txt "`testname' of model assumptions:" ///
										_col(`=`W'-35') "Chi2(`df_`test'') test statistic" ///
										_col(`=`W'-10') as res %10.4f `chi2_`test''
									di as txt _col(`=`W'-35') as txt "p-value" ///
										_col(`=`W'-10') as res %10.4f `p_`test''
								}
								capture confirm number `_shaus_seU'
								if !_rc & strpos(" `_tests_done' "," hausman ") & !missing(`_shaus_seU') {
									di as txt _col(`=`W'-35') "(unrestricted elast." ///
										_col(`=`W'-10') as res %10.4f `_shaus_eU'
									di as txt _col(`=`W'-35') " std. err.)" ///
										_col(`=`W'-10') as res %10.4f `_shaus_seU'
								}
								di as txt "{hline `W'}"
							}
							if `docontrast' {
								tempname _cbb
								matrix `_cbb' = e(b)
								local _cwstub = strlen("`e(depvar)'")
								local _ccn : colnames `_cbb'
								local _ceq : coleq `_cbb'
								foreach x of local _ccn {
									local _cwstub = max(`_cwstub', strlen("`x'"))
								}
								foreach x of local _ceq {
									local _cwstub = max(`_cwstub', strlen("`x'"))
								}
								local _cwstub = max(`_cwstub', 12)
								local _cW = `_cwstub' + 67
								local _cestlab = cond(`estimator'==1,"estimator 1", ///
									cond(`estimator'==2,"estimator 2", ///
									"estimator 3 (`=cond(`useconstant',"constant","exact")', `=cond(`nosplit',"poolmass","splitmass")')"))
								if `_c_ok' {
									di as txt "Elasticity contrast" ///
										_col(`=`_cW'-42') "(reported: `_cestlab';  reference: estimator 3, exact, splitmass)"
									di as txt _col(`=`_cW'-35') "reported elasticity" ///
										_col(`=`_cW'-10') as res %10.4f `_ce0'
									di as txt _col(`=`_cW'-35') "reference elasticity" ///
										_col(`=`_cW'-10') as res %10.4f `_cref0'
									di as txt _col(`=`_cW'-35') "difference" ///
										_col(`=`_cW'-10') as res %10.4f `_d0'
									di as txt _col(`=`_cW'-35') "paired bootstrap SE" ///
										_col(`=`_cW'-10') as res %10.4f `_c_sed'
									di as txt _col(`=`_cW'-35') "corr(reported, reference)" ///
										_col(`=`_cW'-10') as res %10.4f `_c_corr'
									di as txt _col(`=`_cW'-35') "95% CI for difference" ///
										_col(`=`_cW'-21') as res %9.4f `_c_lo' as txt " ," as res %9.4f `_c_hi'
									if !missing(`_c_z') {
										di as txt _col(`=`_cW'-35') "H0: difference = 0    z" ///
											_col(`=`_cW'-10') as res %10.4f `_c_z'
										di as txt _col(`=`_cW'-35') "p-value (normal / bootstrap)" ///
											_col(`=`_cW'-10') as res %10.4f `_c_pnorm' as txt " /" as res %8.4f `_c_ppct'
									}
									else {
										di as txt _col(`=`_cW'-35') "p-value (bootstrap)" ///
											_col(`=`_cW'-10') as res %10.4f `_c_ppct'
									}
									di as txt "{hline `_cW'}"
									di as txt "`_c_nused'/`bootreps' reps paired.  The difference is this estimator's approximation and"
									di as txt "finite-sample gap from the efficient estimator on this sample -- not a signed bias, and"
									di as txt "not comparable across estimators.  See {help polbunch##contrast:help polbunch}."
								}
								else {
									di as txt "Elasticity contrast (reference: estimator 3, exact, splitmass) -- not computed"
									di as txt "  (`_c_why')"
									di as txt "{hline `_cW'}"
								}
							}
							// Analytical bias of the fitted estimator (suppressed
							// by nobias).  polbunchbias reads the counterfactual
							// polynomial, tax rates, window and elasticity from
							// e(); we pass the response-inversion and mass axes
							// explicitly so they match how this estimate was run.
							// Estimators 0/3 are consistent under exact+splitmass,
							// and estimator 4 (Saez) with constant applies a flat
							// counterfactual in both the mass and inversion steps
							// -- both are zero-bias by construction, polbunchbias
							// errors there, and the capture swallows it silently.
							if "`nobias'" == "" & inlist(`estimator', 0, 1, 2, 3, 4) & ///
									"`transform'" != "notransform" & "`t0'" != "" & "`t1'" != "" {
								local _binvopt = cond(`useconstant',"constant","exact")
								local _bmassopt = cond(`nosplit',"poolmass","splitmass")
								local _biteropt "iterate"
								if "`noiterate'" != "" local _biteropt ""
								capture quietly polbunchbias, `_binvopt' `_bmassopt' `_biteropt'
								local _pbbias_ok = (_rc == 0)
								capture confirm number 1   // clear _rc left by a failed polbunchbias (it exits 198 on the zero-bias axes for estimators 0/3)
								if `_pbbias_ok' {
									// Capture r() before any command can overwrite it
									local _bh  = r(bias_h)
									local _bsl = r(bias_slope)
									local _bla = r(bias_lambda)
									local _bB  = r(bias_B)
									local _br  = r(bias_response)
									local _bs  = r(bias_shift)
									local _be  = r(bias_elasticity)
									local _bK  = r(polynomial)
									local _bKfull = r(bias_polyfull)
									local _bl2 = r(bias_l2err)
									local _bdiv = r(iterate_diverged)
									local _biterK  = r(bias_iter_polynomial)
									local _biterl2 = r(bias_iter_l2err)
									local _biterbe = r(bias_iter_elasticity)
									local _biterbs = r(bias_iter_shift)
									local _biterbB = r(bias_iter_B)
									local _biterbr = r(bias_iter_response)
									local _buniterK  = r(bias_uniter_polynomial)
									local _buniterbe = r(bias_uniter_elasticity)
									local _buniterbs = r(bias_uniter_shift)
									local _buniterbB = r(bias_uniter_B)
									local _buniterbr = r(bias_uniter_response)
									tempname _bm _bbeta
									matrix `_bm'    = r(b)
									matrix `_bbeta' = r(bias_beta)

									// Compute table width matching eret di
									tempname _btab
									matrix `_btab' = e(b)
									local _bstub = strlen("`e(depvar)'")
									local _bcols : colnames `_btab'
									local _beqs  : coleq `_btab'
									foreach _bx of local _bcols {
										local _bstub = max(`_bstub', strlen("`_bx'"))
									}
									foreach _bx of local _beqs {
										local _bstub = max(`_bstub', strlen("`_bx'"))
									}
									local _bstub = max(`_bstub', 18)
									local _bW    = `_bstub' + 67

									// Append to eret di table: title + one row per
									// column of r(b) (b0..bK then the estimands).
									local _bnames : colnames `_bm'
									local _bncol  = colsof(`_bm')
									local _bnm1 : word 1 of `_bnames'

									// no reliable bias even after the order fallback:
									// one honest line instead of a table of dots
									if missing(`_bh') {
										di as txt "Bias of the polynomial bunching estimator:" ///
											_col(`=`_bW'-35') as txt "not identified on this window"
										di as txt "{hline `_bW'}"
									}
									else {
									local _bordtxt "order `_bK'"
									if `_buniterK' < . ///
										local _bordtxt "order `_bK', converged iterate; degree-`_buniterK' plug-in did not converge"
									else if `_bKfull' < . & `_bK' < `_bKfull' ///
										local _bordtxt "order `_bK', fallback from `_bKfull'"
									di as txt "Bias of the polynomial bunching estimator (`_bordtxt'):" ///
										_col(`=`_bW'-35') as txt "`_bnm1'" ///
										_col(`=`_bW'-10') as res %10.6g `_bm'[1,1]
									forvalues _bj = 2/`_bncol' {
										local _bnm : word `_bj' of `_bnames'
										di as txt _col(`=`_bW'-35') as txt "`_bnm'" ///
											_col(`=`_bW'-10') as res %10.6g `_bm'[1,`_bj']
									}
									di as txt "{hline `_bW'}"

									if `_buniterK' < . {
										di as text "Note: iterate did not converge at degree `_buniterK' (plug-in elasticity bias"
										di as text "      " %9.4f `_buniterbe' " there); PROMOTED to the converged, self-consistent"
										di as text "      fit at degree `_bK' shown above.  The plug-in number is kept, not discarded"
										di as text "      -- see e(bias_uniter_*).  Specify {cmd:nopromote} to keep it as primary instead."
									}
									else if `_bKfull' < . & `_bK' < `_bKfull' {
										di as text "Note: the degree-`_bKfull' fitting design is too ill-conditioned on this window;"
										di as text "      the bias above is computed at degree `_bK' (best L2 fit to the degree-`_bKfull'"
										di as text "      counterfactual, relative change " %6.4f `_bl2' ///
											cond(`_bl2' < 0.02, " -- the same quantity).", " -- treat as indicative).")
									}

									if "`_bdiv'" == "1" {
										di as error "Warning: the analytical-bias iterate loop did not converge; the un-iterated"
										di as error "         bias is reported above (e(bias_iterate_diverged)=1)."
									}

									if `_biterK' < . {
										di as text "Note: iterate did not converge at degree `_bK', but DOES converge at degree"
										di as text "      `_biterK' (a genuinely self-consistent fixed point; L2 change " %6.4f `_biterl2' "):"
										di as text "      elasticity bias " %9.4f `_biterbe' " there, vs " %9.4f `_be' ///
											" un-iterated at degree `_bK' above (nopromote: not promoted)."
									}

									di as text "Note: this bias is measured against the counterfactual that"
									di as text "      estimator `estimator' itself assumes; each estimator assumes a"
									di as text "      different counterfactual, so these figures are NOT directly"
									di as text "      comparable across estimators.  For a cross-estimator"
									di as text "      comparison, evaluate every estimator against one common"
									di as text "      counterfactual density."
									}

									// Store in e()
									ereturn scalar bias_iterate_diverged = cond("`_bdiv'"=="", ., real("`_bdiv'"))
									ereturn scalar bias_h           = `_bh'
									ereturn scalar bias_slope       = `_bsl'
									ereturn scalar bias_lambda      = `_bla'
									ereturn scalar bias_B           = `_bB'
									ereturn scalar bias_response    = `_br'
									ereturn scalar bias_shift       = `_bs'
									ereturn scalar bias_elasticity  = `_be'
									ereturn scalar bias_polynomial  = `_bK'
									if `_bKfull' < . ereturn scalar bias_polyfull = `_bKfull'
									if `_bl2' < .    ereturn scalar bias_l2err    = `_bl2'
									if `_biterK' < . {
										ereturn scalar bias_iter_polynomial = `_biterK'
										ereturn scalar bias_iter_l2err      = `_biterl2'
										ereturn scalar bias_iter_elasticity = `_biterbe'
										ereturn scalar bias_iter_shift      = `_biterbs'
										ereturn scalar bias_iter_B          = `_biterbB'
										ereturn scalar bias_iter_response   = `_biterbr'
									}
									if `_buniterK' < . {
										ereturn scalar bias_uniter_polynomial = `_buniterK'
										ereturn scalar bias_uniter_elasticity = `_buniterbe'
										ereturn scalar bias_uniter_shift      = `_buniterbs'
										ereturn scalar bias_uniter_B          = `_buniterbB'
										ereturn scalar bias_uniter_response   = `_buniterbr'
									}
									ereturn matrix bias_beta        = `_bbeta'
									ereturn matrix bias             = `_bm'
								}
							}
}
					
					if "`vce'"!="none" set coeftabresults `coeftabresults'

					/*
						savebins() -- write the histogram and the raw bin
						co-visitation matrix to a dataset in the shape that
						vce(cluster stub) reads back:
						    freq  midpoint  m1 m2 ... mJ
						with the cluster count / cutoff / bw kept as dataset
						characteristics.  Done last: it clears the data, and
						polbunch's own preserve restores the caller's data on
						exit.
					*/
					if `"`savebins'"' != "" {
						local _sbJ = rowsof(`_Mraw')
						tempname _sbfr
						frame create `_sbfr'
						frame `_sbfr' {
							svmat double `_sbhist', names(col)
							svmat double `_Mraw', name(cv)
							order freq midpoint cv*
							char _dta[polbunch_nclusters] `_nclust'
							char _dta[polbunch_cutoff]    `cutoff_orig'
							char _dta[polbunch_bw]        `bw_orig'
							label data "polbunch histogram + bin co-visitation matrix"
							save `"`savebinsfile'"', `savebinsrepl'
						}
						frame drop `_sbfr'
						noi di as txt `"(histogram + `_sbJ' x `_sbJ' co-visitation matrix written to `savebinsfile'.dta; `_nclust' clusters)"'
						noi di as txt `"      read back with:  polbunch freq midpoint, ... vce(cluster cv)"'
					}
					}

				end

		cap program drop _pb_parse_clustersub
		program define _pb_parse_clustersub, sclass
			syntax , [ NCLusters(integer -1) ]
			sreturn clear
			sreturn local nclusters = cond(`nclusters' < 0, ".", "`nclusters'")
		end

		cap program drop _pb_parse_vcesub
		program define _pb_parse_vcesub, sclass
			syntax , [ MULTinomial RESidual WILD BAyesian ///
				NORMal BC PERCentile WILDWeights(string) REPS(integer -1) SEED(string) ]

			sreturn clear
			sreturn local boottype = cond("`residual'"!="","residual", ///
				cond("`wild'"!="","wild", ///
				cond("`bayesian'"!="","bayesian","multinomial")))
			sreturn local bootci = cond("`bc'"!="","bc", ///
				cond("`percentile'"!="","percentile","normal"))
			sreturn local wildwt "`wildweights'"
			sreturn local reps   "`reps'"
			sreturn local seed   "`seed'"
		end



		program define bunch_profile, eclass
			version 16.0

			syntax varlist(min=4 max=4 numeric) [if] [in] , ///
				ESTimator(integer) ///
				K(integer) ///
				CUTOFF_orig(real) ///
				BW_orig(real) ///
				cutoff_est(real) ///
				bw_est(real) ///
				L(integer) ///
				H(integer) ///
				zbar_est(real) ///
				zl_excl_orig(real) ///
				zh_excl_orig(real) ///
				zl_excl_est(real) ///
				zh_excl_est(real) ///
				[ nonormalize LOG vce(string) initdelta(real 0.05) DELTAmax(real 1) positive HCType(real -1) MASScorr(real 1) COVis(name) NCLust(real 0) nosplit(integer 1) ]

			gettoken yvar rest : varlist
			gettoken zvar rest : rest
			gettoken sidevar bunch : rest

			if !inlist(`estimator', 0, 1, 2, 3) {
				noi di as err "estimator() must be 0, 1, 2, or 3"
				exit 198
			}

		   marksample touse, novarlist
		   replace `touse' = 0 if missing(`yvar') | missing(`zvar') | missing(`bunch')

			local normalized0 = ("`normalize'" != "nonormalize")
			local islog0      = ("`log'"        != "")
			local positive0 = ("`positive'" != "")
			local dovar0 = ("`vce'"=="analytic")
			
			tempvar y_t z_t side_t bunch_t

			gen double `y_t'     = `yvar'   if `touse'
			gen double `z_t'     = `zvar'   if `touse'
			gen double `side_t'  = `sidevar' if `touse'
			gen double `bunch_t' = `bunch'  if `touse'

			tempname b V Gstack mustack ystack

			capture scalar drop r_weakid_profile r_nbasin_profile

			mata: profile_run( ///
				"`y_t'", ///
				"`z_t'", ///
				"`side_t'", ///
				"`bunch_t'", ///
				`cutoff_orig', ///
				`bw_orig', ///
				`cutoff_est', ///
				`bw_est', ///
				`k', ///
				`estimator', ///
				`islog0', ///
				`zl_excl_orig', ///
				`zh_excl_orig', ///
				`zl_excl_est', ///
				`zh_excl_est', ///
				`zbar_est', ///
				`dovar0', ///
				`initdelta', ///
				`positive0', ///
				`hctype', ///
				`masscorr', ///
				`deltamax', ///
				"`covis'", ///
				`nclust', ///
				`nosplit' ///
			)

			matrix `b' = r_b_profile
			matrix `V' = r_V_profile
			local phihat = .
			capture confirm scalar r_phi_profile
			if !_rc local phihat = r_phi_profile

			/* was the robust/hc2/hc3 sandwich non-PSD and clipped? (see
			   variance_robust) -- absent (missing) under vce(analytic)
			   or vce(cluster), which never set this scalar. */
			local psdclip = .
			capture confirm scalar r_robust_psdclip
			if !_rc local psdclip = r_robust_psdclip

			/* estimator 2/3 profile: weak-identification flag + basin count */
			local weakid0 = 0
			capture confirm scalar r_weakid_profile
			if !_rc local weakid0 = r_weakid_profile
			local nbasin0 = .
			capture confirm scalar r_nbasin_profile
			if !_rc local nbasin0 = r_nbasin_profile

			local h0names
			forvalues j = 1/`k' {
				local term "c.`zvar'"
				if `j' > 1 {
					forvalues r = 2/`j' {
						local term "`term'#c.`zvar'"
					}
				}
				local h0names `h0names' `term'
			}
			local h0names `h0names' _cons

			local cnames
			local eqnames

			if `estimator' == 0 {
				local cnames `h0names' `h0names' B
				forvalues j = 1/`=`k'+1' {
					local eqnames `eqnames' h0
				}
				forvalues j = 1/`=`k'+1' {
					local eqnames `eqnames' h1
				}
				local eqnames `eqnames' bunching
			}
			else if `estimator' == 1 {
				local cnames `h0names' B
				forvalues j = 1/`=`k'+1' {
					local eqnames `eqnames' h0
				}
				local eqnames `eqnames' bunching
			}
			else {
				local cnames `h0names' delta
				forvalues j = 1/`=`k'+1' {
					local eqnames `eqnames' h0
				}
				local eqnames `eqnames' bunching
			}

			matrix colnames `b' = `cnames'
			matrix coleq    `b' = `eqnames'

			if `dovar0' {
				matrix rownames `V' = `cnames'
				matrix colnames `V' = `cnames'
				matrix roweq    `V' = `eqnames'
				matrix coleq    `V' = `eqnames'

				ereturn post `b' `V', esample(`touse')
			}
			else {
				ereturn post `b', esample(`touse')
			}

			ereturn local cmd "bunch_profile"
			ereturn local depvar "`yvar'"
			ereturn local zvar "`zvar'"

			ereturn scalar estimator   = `estimator'
			if `phihat'<. ereturn scalar dispersion = `phihat'
			if `psdclip'<. ereturn scalar vce_psdclip = `psdclip'
			if inlist(`estimator',2,3) {
				ereturn scalar delta_weakid = `weakid0'
				if `nbasin0'<. ereturn scalar delta_nbasin = `nbasin0'
			}
			ereturn scalar K           = `k'
			ereturn scalar cutoff_orig = `cutoff_orig'
			ereturn scalar bw_orig     = `bw_orig'
			ereturn scalar L           = `l'
			ereturn scalar H           = `h'
			ereturn scalar normalized  = `normalized0'
			ereturn scalar islog       = `islog0'
			ereturn scalar positive = `positive0'
			

			if `dovar0' {
				ereturn local vcetype "Analytic"
				ereturn local properties "b V"

				matrix `Gstack'  = r_G_stack
				matrix `mustack' = r_mu_stack
				matrix `ystack'  = r_ystack

				ereturn matrix G_stack  = `Gstack'
				ereturn matrix mu_stack = `mustack'
				ereturn matrix y_stack  = `ystack'

				_pb_gof_post
			}
		end

		/* posts the reference-region goodness-of-fit scalars from the
		   Stata objects r_gof / r_gof_rmse / r_gof_massresid last set by
		   profile_run() or saez_run().  eclass, no -ereturn post-, so it
		   adds e(*) scalars to the estimation results already in place. */
		cap program drop _pb_gof_post
		program define _pb_gof_post, eclass
			capture confirm matrix r_gof
			if _rc exit
			tempname G
			matrix `G' = r_gof
			if colsof(`G') < 13 exit
			ereturn scalar gof_nbins     = `G'[1,1]
			ereturn scalar gof_np        = `G'[1,2]
			ereturn scalar gof_df        = `G'[1,3]
			ereturn scalar deviance      = `G'[1,4]
			ereturn scalar pearson_x2    = `G'[1,5]
			ereturn scalar gof_ll        = `G'[1,6]
			ereturn scalar deviance_null = `G'[1,7]
			ereturn scalar gof_ll_null   = `G'[1,8]
			ereturn scalar r2_dev        = `G'[1,9]
			ereturn scalar aic           = `G'[1,10]
			ereturn scalar bic           = `G'[1,11]
			ereturn scalar qaic          = `G'[1,12]
			ereturn scalar qaicc         = `G'[1,13]
			if `G'[1,3] > 0 & `G'[1,3] < . {
				ereturn scalar deviance_p   = chi2tail(`G'[1,3], `G'[1,4])
				ereturn scalar pearson_x2_p = chi2tail(`G'[1,3], `G'[1,5])
			}
			capture confirm scalar r_gof_rmse
			if !_rc ereturn scalar gof_rmse = r_gof_rmse
			capture confirm scalar r_gof_massresid
			if !_rc ereturn scalar gof_massresid = r_gof_massresid

			foreach _s in below above {
				capture confirm matrix r_gof_`_s'
				if _rc continue
				tempname S
				matrix `S' = r_gof_`_s'
				if colsof(`S') < 5 continue
				if missing(`S'[1,4]) continue
				ereturn scalar gof_df_`_s'     = `S'[1,3]
				ereturn scalar deviance_`_s'   = `S'[1,4]
				ereturn scalar pearson_x2_`_s' = `S'[1,5]
				if `S'[1,3] > 0 & `S'[1,3] < . {
					ereturn scalar dispersion_`_s' = `S'[1,5] / `S'[1,3]
					ereturn scalar deviance_`_s'_p = chi2tail(`S'[1,3], `S'[1,4])
				}
			}
		end


		program define bunch_transform, eclass
			version 16.0

			syntax varname, ///
				ESTimator(integer) ///
				K(integer) ///
				CUTOFFORIG(real) ///
				CUTOFFEST(real) ///
				BWORIG(real) ///
				BWEST(real) ///
				XSCALE(real) ///
				ZLEXCL(real) ///
				ZHEXCL(real) ///
				ZLEXCLEST(real) ///
				ZHEXCLEST(real) ///
				LOW(integer) ///
				HIGH(integer) ///
				[ LOG CONSTANT T0(numlist min=1 max=1) T1(numlist min=1 max=1) nograd nonormalize ZBAR(real 0) MASSOBS(real 0) NOSPLIT(integer 0) WEAKid(integer 0) ]

				loc z `varlist'
				
			if "`t0'" != "" {
				local t0 : word 1 of `t0'
			}
			if "`t1'" != "" {
				local t1 : word 1 of `t1'
			}

			/* Existing e(b) is required */
			capture confirm matrix e(b)
			if _rc {
				di as err "e(b) not found"
				exit 301
			}

			tempname theta Vtheta bnew Gnew Vnew

			matrix `theta' = e(b)

			/* Use delta-method VCE only if e(V) exists and nograd is not specified */
			local dograd = 0
			if "`nograd'" == "" {
				capture confirm matrix e(V)
				if !_rc {
					matrix `Vtheta' = e(V)
					local dograd = 1
				}
			}

			/* Flags */
			local islog       = ("`log'"      != "")
			local constant0   = ("`constant'" != "")
			/*
				nosplit0: 1 = poolmass  (B = M - int_{zL}^{zH} h0)
				          0 = splitmass (B = M - int_{zL}^{z*} h0
				                              - int_{z*}^{zH} h1)
				Ignored for estimator 1 (h0 == h1).
			*/
			local nosplit0 = `nosplit'
			/*
				A structural delta column is reported for estimator 2
				always, and for estimator 3 whenever the reported response
				is NOT the pure structural delta -- i.e. under constant
				(constant-density inversion) or under poolmass (response
				backed out from the pooled reduced-form B).
			*/
			local e3delta = (`estimator' == 3 & (`constant0' | `nosplit0'))
			local hastax0 = ("`t0'" != "" & "`t1'" != "")
			if `hastax0' {
				if "`t0'" == "" | "`t1'" == "" {
					di as err "options t0() and t1() are required when tax options are used"
					exit 198
				}
			}
			else {
				if "`t0'" == "" local t0 = 0
				if "`t1'" == "" local t1 = 0
			}

			/* Preserve e(sample), if present */
			tempvar touse
			capture gen byte `touse' = e(sample)
			local has_esample = !_rc

			mata: bunch_transform( ///
				st_matrix("`theta'"), ///
				`estimator', ///
				`k', ///
				`cutofforig', ///
				`cutoffest', ///
				`bworig', ///
				`bwest', ///
				`xscale', ///
				`islog', ///
				`constant0', ///
				`hastax0', ///
				`t0', ///
				`t1', ///
				`dograd', ///
				`zbar', ///
				`massobs', ///
				`nosplit0', ///
				`zlexcl', ///
				`zhexcl', ///
				`zlexclest', ///
				`zhexclest', ///
				`weakid' ///
			)

			matrix `bnew' = b_bunchcalc
			local hasresp = 1

			capture confirm scalar b_bunchcalc_hasresp
			if !_rc {
				local hasresp = scalar(b_bunchcalc_hasresp)
			}

			if `dograd' {
				matrix `Gnew' = G_bunchcalc

				if colsof(`Gnew') != colsof(`theta') {
					noi di as err "conformability error: colsof(G) != colsof(e(b))"
					exit 503
				}

				if rowsof(`Vtheta') != colsof(`theta') | colsof(`Vtheta') != colsof(`theta') {
					di as err "conformability error: e(V) is not compatible with e(b)"
					exit 503
				}

				matrix `Vnew' = `Gnew' * `Vtheta' * `Gnew''
			}

		   /* Coefficient names with equations */
			local h0names
			local h1names

			forvalues j = 1/`k' {
				if `j' == 1 {
					local term c.`z'
				}
				else {
					local term `term'#c.`z'
				}

				local h0names `h0names' `term'
				local h1names `h1names' `term'
			}

			local h0names `h0names' _cons
			local h1names `h1names' _cons

			local cnames ///
				`h0names' ///
				`h1names' ///
				number_bunchers ///
				excess_mass

			if (`estimator' == 2 | `e3delta') {
				local cnames `cnames' delta
			}

			if `hasresp' {
				local cnames `cnames' ///
					shift ///
					marginal_response

				if `hastax0' {
					local cnames `cnames' elasticity
				}
			}
			else if `weakid' {
				noi di as text "Note: the response length (delta) for estimator `estimator' is weakly identified for this"
				noi di as text "      window / polynomial order -- the profile search hit the deltamax() bound or"
				noi di as text "      found competing minima. The counterfactual, number of bunchers and excess mass"
				noi di as text "      are reported; shift, marginal response and elasticity are withheld. Check the"
				noi di as text "      minimum-distance / Hausman test and try a different window or polynomial()."
			}
			else {
				noi di as text "Note: Could not find real root to solve the polynomial. Consider using the constant approximation."
			}

			if wordcount("`cnames'") != colsof(`bnew') {
				di as err "internal error: coefficient names do not match transformed b"
				di as err "number of names = " wordcount("`cnames'")
				di as err "colsof(b)       = " colsof(`bnew')
				exit 503
			}

			matrix colnames `bnew' = `cnames'


			/* Equation names */
			local eqnames

			forvalues j = 1/`=`k'+1' {
				local eqnames `eqnames' h0
			}

			forvalues j = 1/`=`k'+1' {
				local eqnames `eqnames' h1
			}

			local eqnames `eqnames' bunching bunching

			if (`estimator' == 2 | `e3delta') {
				local eqnames `eqnames' bunching
			}

			if `hasresp' {
				local eqnames `eqnames' bunching bunching

				if `hastax0' {
					local eqnames `eqnames' bunching
				}
			}

			matrix coleq `bnew' = `eqnames'

			if `dograd' {
				matrix rownames `Vnew' = `cnames'
				matrix colnames `Vnew' = `cnames'
				matrix roweq    `Vnew' = `eqnames'
				matrix coleq    `Vnew' = `eqnames'

				matrix rownames `Gnew' = `cnames'
				matrix roweq    `Gnew' = `eqnames'
			}

			/* Post transformed results */
			if `dograd' {
				if `has_esample' {
					ereturn post `bnew' `Vnew', esample(`touse')
				}
				else {
					ereturn post `bnew' `Vnew'
				}

				ereturn matrix G = `Gnew'
				ereturn local vcetype "delta method"
				ereturn local properties "b V"
			}
			else {
				if `has_esample' {
					ereturn post `bnew', esample(`touse')
				}
				else {
					ereturn post `bnew'
				}
			}

			ereturn local cmd "bunch_transform"
			/*
				hasresp: whether shift/marginal_response/elasticity are
				identified for this draw (the exact Naive/Chetty/Saez
				quadratic had a real root). When 0, those three
				coefficients are absent from e(b) entirely and
				"_b[bunching:elasticity]" (etc.) errors with "not found"
				rather than returning a missing value -- Stata's
				"ereturn post"/"ereturn repost" both unconditionally
				refuse any b vector containing missing entries, so the
				width has to vary rather than posting them as missing.
				Check e(hasresp) BEFORE touching those coefficients
				instead of wrapping the access in capture.
			*/
			ereturn scalar hasresp = `hasresp'

			ereturn scalar estimator   = `estimator'
			ereturn scalar K           = `k'
			ereturn scalar cutoff_orig = `cutofforig'
			ereturn scalar cutoff_est  = `cutoffest'
			ereturn scalar bw_orig     = `bworig'
			ereturn scalar bw_est      = `bwest'
			ereturn scalar xscale      = `xscale'
			ereturn scalar L           = `low'
			ereturn scalar H           = `high'
			ereturn scalar islog       = `islog'
			ereturn scalar constant    = `constant0'
			ereturn scalar nosplit     = `nosplit0'
			ereturn scalar hastax      = `hastax0'

			if `hastax0' {
				ereturn scalar t0 = `t0'
				ereturn scalar t1 = `t1'
			}
		end
		

capture program drop saez_transform
program define saez_transform, eclass
    version 16.0

    syntax , ZSTAROrig(numlist max=1) BWOrig(numlist max=1) ///
        [T0(numlist max=1) T1(numlist max=1) log nograd constant]

    tempname b0 V0 theta Vtheta b G V

    matrix `b0' = e(b)
    if (colsof(`b0') < 3) {
        di as err "e(b) must contain theta=(h0:_cons,h1:_cons,B) in columns 1..3"
        exit 503
    }
	
    local islog       = ("`log'" != "")
    /*
        syntax turns "nograd" into a toggle option whose macro is `grad'
        (empty, or "nograd"), so `nograd' itself is never set.  Only take
        the delta-method path when the caller did NOT pass nograd AND an
        unrestricted e(V) is actually present (mirrors bunch_transform).
    */
    local dograd = ("`grad'`nograd'" != "nograd")
    if (`dograd') {
        capture confirm matrix e(V)
        if (_rc) local dograd = 0
    }
    local useconstant = ("`constant'" != "")

    local zstarorig : word 1 of `zstarorig'
    local bworig    : word 1 of `bworig'

    if (`zstarorig' <= 0) {
        di as err "zstarorig() must be positive"
        exit 198
    }
    if (`bworig' <= 0) {
        di as err "bworig() must be positive"
        exit 198
    }

    local hastax = 0
    if ("`t0'" != "" | "`t1'" != "") {
        if ("`t0'" == "" | "`t1'" == "") {
            di as err "t0() and t1() must be specified together"
            exit 198
        }
        local t0 : word 1 of `t0'
        local t1 : word 1 of `t1'
        local hastax = 1
    }
    else {
        local t0 = .
        local t1 = .
    }

    matrix `theta' = `b0'[1,1..3]
    matrix colnames `theta' = h0:_cons h1:_cons bunching:number_bunchers

    if (`dograd') {
        capture matrix `Vtheta' = e(V)
        if (_rc) {
            di as err "e(V) not found; specify nograd"
            exit 111
        }
    }

    mata: st_matrix("`b'", saez_transform( ///
        st_matrix("`theta'"), ///
        `zstarorig', `bworig', `t0', `t1', ///
        `islog', `hastax', `dograd', `useconstant', "`G'" ///
    ))

    local outnames h0:_cons h1:_cons bunching:number_bunchers ///
        bunching:excess_mass bunching:shift bunching:marginal_response
    if (`hastax') local outnames `outnames' bunching:elasticity

    matrix colnames `b' = `outnames'

    if (`dograd') {
        matrix rownames `G' = `outnames'
        matrix colnames `G' = h0:_cons h1:_cons bunching:number_bunchers

        matrix `V' = `G' * `Vtheta' * `G''
        matrix rownames `V' = `outnames'
        matrix colnames `V' = `outnames'

        ereturn post `b' `V'
        ereturn matrix G = `G'
        ereturn matrix Vtheta = `Vtheta'
    }
    else {
        ereturn post `b'
    }

    ereturn matrix theta = `theta'
    ereturn scalar zstarorig = `zstarorig'
    ereturn scalar bworig    = `bworig'
    ereturn scalar islog     = `islog'
    ereturn scalar hastax    = `hastax'

    if (`hastax') {
        ereturn scalar t0 = `t0'
        ereturn scalar t1 = `t1'
    }

    ereturn local cmd "saez_transform"
	ereturn local properties "b V"
    ereturn display
end

cap program drop bunch_saez
	program define bunch_saez, eclass
		version 16.0

		syntax varlist(min=4 max=4 numeric) [if] [in], ///
			CUTOFF_orig(real) ///
			BW_orig(real) ///
			ZL_excl_orig(real) ///
			ZH_excl_orig(real) ///
			[vce(string) HCType(real -1) MASScorr(real 1) COVis(name) NCLust(real 0)]

		gettoken yvar rest : varlist
		gettoken zvar rest : rest
		gettoken sidevar rest : rest
		gettoken bunchvar : rest

		marksample touse, novarlist
		replace `touse' = 0 if missing(`yvar') | missing(`zvar') | missing(`bunchvar')

		
		quietly count if `touse' & `bunchvar' == 0 & missing(`sidevar')
		if r(N) > 0 {
			di as err "Some non-excluded bins cannot be classified as left or right of cutoff."
			exit 498
		}

		quietly count if `touse' & `bunchvar' == 0 & `sidevar' == -1
		if r(N) == 0 {
			di as err "no left reference bins found for Saez estimator"
			exit 498
		}

		quietly count if `touse' & `bunchvar' == 0 & `sidevar' == 1
		if r(N) == 0 {
			di as err "no right reference bins found for Saez estimator"
			exit 498
		}

		quietly count if `touse' & `bunchvar' > 0
		if r(N) == 0 {
			di as err "no excluded/bunching bins found for Saez estimator"
			exit 498
		}

		local width_excl = r(N)

		/*
			Saez counterfactual weights inside the excluded region.

			a0 is the width, in bins, of the excluded interval below the cutoff.
			a1 is the width, in bins, of the excluded interval above the cutoff.

			This gives:
				Hstar = a0*h0 + a1*h1 + B
			hence:
				B = Hstar - a0*h0 - a1*h1

			For symmetric excluded regions this collapses to:
				a0 = a1 = 0.5*width_excl.
		*/
		local a0 = (`cutoff_orig' - `zl_excl_orig') / `bw_orig'
		local a1 = (`zh_excl_orig' - `cutoff_orig') / `bw_orig'

		if `a0' < -1e-8 | `a1' < -1e-8 {
			di as err "invalid Saez excluded-region weights"
			di as err "a0 = " %12.8f `a0' ", a1 = " %12.8f `a1'
			exit 498
		}

		if abs((`a0' + `a1') - `width_excl') > 1e-6 {
			di as err "internal error: Saez weights do not sum to excluded-region width"
			di as err "a0 + a1 = " %12.8f (`a0' + `a1') ///
				", excluded bins = " %12.8f `width_excl'
			exit 498
		}

		/*
			Clean tiny floating-point artifacts.
		*/
		if abs(`a0') < 1e-10 local a0 = 0
		if abs(`a1') < 1e-10 local a1 = 0

		local dovar0 = ("`vce'"=="analytic")
		
		tempvar y_t side_t bunch_t
		gen double `y_t'     = `yvar'     if `touse'
		gen double `side_t'  = `sidevar'  if `touse'
		gen double `bunch_t' = `bunchvar' if `touse'
		

		mata: saez_run("`y_t'", "`side_t'", "`bunch_t'", `a0', `a1', `dovar0', `hctype', `masscorr', "`covis'", `nclust')
		local phihat = .
		capture confirm scalar r_phi_saez
		if !_rc local phihat = r_phi_saez

		tempname b V Gstack mustack ystack

		matrix `b' = r_b_saez
		matrix colnames `b' = _cons _cons B
		matrix coleq    `b' = h0 h1 bunching
		
		matrix `V' = r_V_saez

		if `dovar0' {
			matrix rownames `V' = _cons _cons B
			matrix colnames `V' = _cons _cons B
			matrix roweq    `V' = h0 h1 bunching
			matrix coleq    `V' = h0 h1 bunching
			ereturn post `b' `V', esample(`touse')
		}
		else {
			ereturn post `b', esample(`touse')
		}

		if `dovar0' {
			ereturn local vcetype "Analytic"
			ereturn local properties "b V"
			
			matrix `Gstack' = r_G_saez
			matrix `mustack' = r_mu_saez
			matrix `ystack' = r_ystack_saez

			ereturn matrix G_stack = `Gstack'
			ereturn matrix mu_stack = `mustack'
			ereturn matrix y_stack = `ystack'

			_pb_gof_post
		}

		local psdclip = .
		capture confirm scalar r_robust_psdclip
		if !_rc local psdclip = r_robust_psdclip

		ereturn local cmd "bunch_saez"
		ereturn scalar estimator = 4
		if `phihat'<. ereturn scalar dispersion = `phihat'
		if `psdclip'<. ereturn scalar vce_psdclip = `psdclip'
		ereturn scalar cutoff_orig = `cutoff_orig'
		ereturn scalar bw_orig = `bw_orig'
		ereturn scalar zL_excl_orig = `zl_excl_orig'
		ereturn scalar zH_excl_orig = `zh_excl_orig'
		ereturn scalar saez_a0 = `a0'
		ereturn scalar saez_a1 = `a1'
		ereturn scalar saez_width_excl = `width_excl'
	end

	cap program drop polbunch_shausmantest
	program define polbunch_shausmantest, rclass
		version 16.0

		/*
			Scalar Hausman test on the elasticity.  The caller captures the
			inputs from the restricted fit + its bunch_transform and from the
			unrestricted (estimator 0) fit + its bunch_transform:

			  gbcr / gbcu : bunch_transform Jacobian e(G) of the restricted /
			                unrestricted model (elasticity = last row)
			  gr   / gu   : stacked-fit gradient e(G_stack) of each model
			  er   / eu   : the two elasticity point estimates (may be .)
			  ystack      : e(y_stack) -- the shared bin counts
			  muu         : e(mu_stack) of the unrestricted fit (HC meat only)
		*/
		syntax , GBCR(name) GR(name) ER(string) ///
		         GBCU(name) GU(name) EU(string) ///
		         YSTACK(name) ///
		         [ MUU(name) HCType(real -1) COVis(name) NCLust(real 0) ]

		local erv = real("`er'")
		local euv = real("`eu'")

		local muname "."
		if "`muu'" != "" {
			capture confirm matrix `muu'
			if !_rc local muname "`muu'"
		}
		local covname "."
		if "`covis'" != "" {
			capture confirm matrix `covis'
			if !_rc local covname "`covis'"
		}

		capture noisily mata: polbunch_shausman_mata( ///
			"`gbcu'", "`gu'", `euv', ///
			"`gbcr'", "`gr'", `erv', ///
			"`ystack'", ///
			"`muname'", `hctype', "`covname'", `nclust' )

		if _rc {
			return scalar chi2 = .
			return scalar p    = .
			return scalar df   = .
			return scalar seU  = .
			return scalar failcode = _rc
			exit
		}

		return scalar chi2     = r(pb_sh_chi2)
		return scalar p        = r(pb_sh_p)
		return scalar df       = r(pb_sh_df)
		return scalar seU      = r(pb_sh_seU)
		return scalar failcode = r(pb_sh_failcode)
	end



		cap prog drop polbunch_waldtest
		program define polbunch_waldtest, rclass
			version 16.0

			syntax , ///
				ESTimator(integer) ///
				K(integer) ///
				CUTOFFORIG(real) ///
				CUTOFFEST(real) ///
				BWORIG(real) ///
				BWEST(real) ///
				ZBAR(real) ///
				[ nonormalize LOG ]

			if !inlist(`estimator', 1, 2, 3) {
				di as err "polbunch_waldtest only handles estimator(1), estimator(2), or estimator(3)"
				exit 198
			}

			capture confirm matrix e(b)
			if _rc {
				di as err "e(b) not found; post unrestricted estimator(0) before calling polbunch_waldtest"
				exit 301
			}

			capture confirm matrix e(V)
			if _rc {
				di as err "e(V) not found; model-restriction test requires unrestricted VCE"
				exit 301
			}

			tempname b V
			matrix `b' = e(b)
			matrix `V' = e(V)

			local islog0 = ("`log'" != "")

			mata: polbunch_wald_from_unrestricted( ///
				"`b'", ///
				"`V'", ///
				`estimator', ///
				`cutofforig', ///
				`bworig', ///
				`cutoffest', ///
				`bwest', ///
				`k', ///
				`islog0', ///
				`zbar' ///
			)

			tempname chi2 p df deltaU failcode

			scalar `chi2'    = r(pb_wald)
			scalar `p'       = r(pb_p)
			scalar `df'      = r(pb_df)
			scalar `deltaU'  = r(pb_delta_U)
			scalar `failcode' = r(pb_failcode)
			return scalar chi2    = `chi2'
			return scalar p       = `p'
			return scalar df      = `df'
			return scalar delta_U = `deltaU'
			return scalar failcode = `failcode'

		end
		
		cap program drop polbunch_minimumdistancetest
	program define polbunch_minimumdistancetest, rclass
		version 16.0

		syntax , ///
			ESTimator(integer) ///
			K(integer) ///
			CUTOFFORIG(real) ///
			CUTOFFEST(real) ///
			BWORIG(real) ///
			BWEST(real) ///
			ZBAR(real) ///
			[ NONORMALIZE LOG POSitive INITDELTA(real 0.05) ]

		if !inlist(`estimator', 1, 2, 3) {
			di as err "polbunch_minimumdistancetest only handles estimator(1), estimator(2), or estimator(3)"
			return scalar failcode = 198
			return scalar chi2 = .
			return scalar p = .
			return scalar df = .
			return scalar delta = .
			exit 198
		}

		capture confirm matrix e(b)
		if _rc {
			di as err "e(b) not found; post unrestricted estimator(0) before calling polbunch_minimumdistancetest"
			return scalar failcode = 301
			return scalar chi2 = .
			return scalar p = .
			return scalar df = .
			return scalar delta = .
			exit 301
		}

		capture confirm matrix e(V)
		if _rc {
			di as err "e(V) not found; minimum-distance test requires unrestricted VCE"
			return scalar failcode = 302
			return scalar chi2 = .
			return scalar p = .
			return scalar df = .
			return scalar delta = .
			exit 301
		}

		tempname b V
		matrix `b' = e(b)
		matrix `V' = e(V)

		local islog0   = ("`log'" != "")
		local positive0 = ("`positive'" != "")

		capture noisily mata: polbunch_mdt_mata( ///
			"`b'", ///
			"`V'", ///
			`estimator', ///
			`cutofforig', ///
			`bworig', ///
			`cutoffest', ///
			`bwest', ///
			`k', ///
			`islog0', ///
			`zbar', ///
			`positive0', ///
			`initdelta' ///
		)

		if _rc {
			return scalar chi2 = .
			return scalar p = .
			return scalar df = .
			return scalar delta = .
			return scalar failcode = _rc
			exit
		}

		tempname chi2 p df delta failcode

		scalar `chi2'    = r(pb_md)
		scalar `p'       = r(pb_md_p)
		scalar `df'      = r(pb_md_df)
		scalar `delta'   = r(pb_md_delta)
		scalar `failcode' = r(pb_md_failcode)

		if missing(`chi2') | `failcode' {
			di as err "Could not compute minimum-distance statistic."
			di as err "failcode = " `failcode'
			return scalar chi2 = .
			return scalar p = .
			return scalar df = .
			return scalar delta = .
			return scalar failcode = `failcode'
			exit 498
		}

		return scalar chi2 = `chi2'
		return scalar p = `p'
		return scalar df = `df'
		return scalar delta = `delta'
		return scalar failcode = 0
	end

		mata:


		// MAIN STRUCTS
		struct hcoef_out {
			real rowvector gamma
			real matrix dgamma_dbeta
			real colvector dgamma_ddelta
		}

		struct hdesign_out {
			real matrix X
			real matrix dXddelta
		}

		struct stack23_out {
			real colvector ystack
			real matrix X
			real matrix G
			real colvector mu
		}

		struct design_out {
			real matrix X              // stacked design wrt linear parameters
			real matrix dXddelta       // derivative of X wrt delta for estimators 2/3
		}


		// -----------------------------------------------------------------------------
		// Basic polynomial helpers
		// -----------------------------------------------------------------------------

		real matrix pbasis(real colvector z, real scalar K)
		{
			real scalar j
			real matrix X

			X = J(rows(z), K+1, 1)

			for (j=1; j<=K; j++) {
				X[,j] = z:^j
			}

			// constant last
			X[,K+1] = J(rows(z), 1, 1)

			return(X)
		}

		real rowvector pbasis_row(real scalar z, real scalar K)
		{
			real scalar j
			real rowvector x

			x = J(1, K+1, 1)

			for (j=1; j<=K; j++) {
				x[j] = z^j
			}

			// constant last
			x[K+1] = 1

			return(x)
		}

		real rowvector intbasis(real scalar a, real scalar b, real scalar K)
		{
			real scalar j
			real rowvector r

			r = J(1, K+1, .)

			for (j=1; j<=K; j++) {
				r[j] = (b^(j+1) - a^(j+1))/(j+1)
			}

			// constant last
			r[K+1] = b - a

			return(r)
		}

	real scalar response_length(
		real scalar delta,
		real scalar cutoff_orig,
		real scalar bw_orig,
		real scalar bw_est,
		real scalar islog
	)
	{
		if (1 + delta <= 0) _error(3498, "delta must be greater than -1")

		if (islog == 1) {
			return(ln(1 + delta) * bw_est / bw_orig)
		}

		return(delta * cutoff_orig * bw_est / bw_orig)
	}
		
	real scalar d_response_length_ddelta(
		real scalar delta,
		real scalar cutoff_orig,
		real scalar bw_orig,
		real scalar bw_est,
		real scalar islog
	)
	{
		if (1 + delta <= 0) _error(3498, "delta must be greater than -1")

		if (islog == 1) {
			return((bw_est / bw_orig) / (1 + delta))
		}

		return(cutoff_orig * bw_est / bw_orig)
	}

	void saez_run(
		string scalar yvar,
		string scalar sidevar,
		string scalar bunchvar,
		real scalar a0,
		real scalar a1,
		real scalar dovar,
		real scalar hctype,
		real scalar masscorr,
		string scalar covisname,
		real scalar nclust
	)
	{
		real colvector y, side, bunch
		real colvector yL, yR, ystack, mu
		real scalar nL, nR, Hstar_obs
		real matrix X, Vout
		real rowvector theta

		y     = st_data(., yvar)
		side  = st_data(., sidevar)
		bunch = st_data(., bunchvar)

		yL = select(y, (bunch :== 0) :& (side :== -1))
		yR = select(y, (bunch :== 0) :& (side :==  1))

		nL = rows(yL)
		nR = rows(yR)

		Hstar_obs = sum(select(y, bunch :> 0))

		ystack = yL \ yR \ Hstar_obs

		X =
			(J(nL, 1, 1), J(nL, 1, 0), J(nL, 1, 0)) \
			(J(nR, 1, 0), J(nR, 1, 1), J(nR, 1, 0)) \
			(a0,           a1,           1)

		theta = qrsolve(X, ystack)'

		mu = X * theta'

		if (dovar == 1) {
			if (covisname != "") {
				Vout = variance_cluster(X, ystack, st_matrix(covisname), nclust)
			}
			else if (hctype >= 0) {
				Vout = variance_robust(X, ystack, mu, hctype, 1, masscorr)
			}
			else {
				Vout = variance_multinomial(X, ystack, 0)
			}
			st_numscalar("r_phi_saez", pearson_phi(ystack[1::(nL+nR)], mu[1::(nL+nR)], cols(X)))
			st_matrix("r_gof", pbx_gof(ystack[1::(nL+nR)], mu[1::(nL+nR)], cols(X)))
			st_numscalar("r_gof_rmse", sqrt(mean((ystack[1::(nL+nR)] :- mu[1::(nL+nR)]):^2)))
			st_numscalar("r_gof_massresid",
				(ystack[nL + nR + 1] > 0
				 ? (ystack[nL + nR + 1] - mu[nL + nR + 1]) / ystack[nL + nR + 1]
				 : .))
			/* reference fit split at the cutoff (see profile_run); Saez
			   fits one level per side, so p = 1 each */
			st_matrix("r_gof_below",
				(nL >= 2 ? pbx_gof(ystack[1::nL], mu[1::nL], 1) : J(1, 13, .)))
			st_matrix("r_gof_above",
				(nR >= 2 ? pbx_gof(ystack[(nL+1)::(nL+nR)], mu[(nL+1)::(nL+nR)], 1) : J(1, 13, .)))
		}
		else {
			Vout = J(3, 3, .)
		}

			st_matrix("r_b_saez", theta)
			st_matrix("r_V_saez", Vout)

		if (dovar == 1) {
			st_matrix("r_G_saez", X)
			st_matrix("r_mu_saez", mu)
			st_matrix("r_ystack_saez", ystack)

			/*
				Full-length fitted-mean vector for the residual bootstrap --
				the Saez counterfactual (piecewise level a0/a1 on each side)
				on the reference bins, missing on the excluded bins.  See the
				matching block in profile_run().
			*/
			{
				real colvector rb_idx, rb_mu

				rb_idx = selectindex((bunch :== 0) :& (side :== -1)) \
				         selectindex((bunch :== 0) :& (side :==  1))
				rb_mu = J(rows(y), 1, .)
				if (rows(rb_idx) == nL + nR & nL + nR >= 1) {
					rb_mu[rb_idx] = mu[1::(nL + nR)]
				}
				if (rows(rb_mu) <= 10000) st_matrix("r_mu_bins", rb_mu)
			}
		}
	}

	real matrix h1_A_matrix(
		real scalar delta,
		real scalar estimator,
		real scalar K,
		real scalar cutoff_orig,
		real scalar bw_orig,
		real scalar cutoff_est,
		real scalar bw_est,
		real scalar islog,
		real scalar deriv
	)
	{
		real scalar Kb, p, j, e
		real scalar scale, s, a, dscale, ds, da
		real scalar apow, apowm1, spow, spowm1, base, dbase
		real matrix A

		if (1 + delta <= 0) {
			_error(3498, "delta must be greater than -1")
		}

		Kb = K + 1
		A  = J(Kb, Kb, 0)

		/*
			Estimator 1: h1 = h0
		*/
		if (estimator == 1) {
			if (deriv) return(J(Kb, Kb, 0))
			return(I(Kb))
		}

		/*
			Estimator 2: h1 = h0 / (1 + delta)
		*/
		if (estimator == 2) {
			if (deriv) return(-I(Kb) / (1 + delta)^2)
			return(I(Kb) / (1 + delta))
		}

		if (estimator != 3) {
			_error(3498, "h1_A_matrix only handles estimators 1, 2, and 3")
		}

		/*
			Estimator 3.

			Level running variable:
				h1(z) = (1+delta) h0(a + (1+delta)z)

			Log running variable:
				h1(z) = h0(a + z)
		*/

		else {
			if (islog == 0) {
				scale  = 1 + delta
				s      = 1 + delta
				dscale = 1
				ds     = 1

				/*
					Original-variable transformation:
						z0 = (1 + delta) * z

					In normalized coordinates x = (z-zmid)/xscale:
						x0 = (1 + delta)*x + delta*zmid/xscale
				*/
				a = delta * (
					cutoff_orig * bw_est / bw_orig
					- cutoff_est
				)

				da = (
					cutoff_orig * bw_est / bw_orig
					- cutoff_est
				)
			}
			else {
				scale  = 1
				s      = 1
				dscale = 0
				ds     = 0

				a  = ln(1 + delta) * bw_est / bw_orig
				da = (bw_est / bw_orig) / (1 + delta)
			}
		}

		/*
			Nonconstant rows. Coefficients are ordered:
				beta_1, ..., beta_K, beta_0
		*/
		for (p = 1; p <= K; p++) {
			for (j = p; j <= K; j++) {
				e = j - p

				if (e == 0) apow = 1
				else        apow = a^e

				spow = s^p
				base = apow * spow

				if (deriv == 0) {
					A[p,j] = scale * comb(j,p) * base
				}
				else {
					dbase = 0

					if (e > 0) {
						if (e == 1) apowm1 = 1
						else        apowm1 = a^(e-1)

						dbase = dbase + e * apowm1 * da * spow
					}

					if (p > 0) {
						if (p == 1) spowm1 = 1
						else        spowm1 = s^(p-1)

						dbase = dbase + apow * p * spowm1 * ds
					}

					A[p,j] = comb(j,p) * (dscale * base + scale * dbase)
				}
			}
		}

		/*
			Constant row.
		*/
		if (deriv == 0) {
			A[Kb,Kb] = scale

			for (j = 1; j <= K; j++) {
				A[Kb,j] = scale * a^j
			}
		}
		else {
			A[Kb,Kb] = dscale

			for (j = 1; j <= K; j++) {
				if (j == 1) apowm1 = 1
				else        apowm1 = a^(j-1)

				A[Kb,j] = dscale * a^j + scale * j * apowm1 * da
			}
		}

		return(A)
	}

		// -----------------------------------------------------------------------------
		// h1 coefficient/design restrictions
		// -----------------------------------------------------------------------------
	struct hcoef_out scalar h1coef_map(
		real rowvector beta,
		real scalar delta,
		real scalar estimator,
		real scalar K,
		real scalar cutoff_orig,
		real scalar bw_orig,
		real scalar cutoff_est,
		real scalar bw_est,
		real scalar islog,
		real scalar dograd
	)
	{
		struct hcoef_out scalar out
		real scalar Kb
		real matrix A, dA

		if (1 + delta <= 0) {
			_error(3498, "delta must be greater than -1")
		}

		Kb = K + 1

		A = h1_A_matrix(
			delta,
			estimator,
			K,
			cutoff_orig,
			bw_orig,
			cutoff_est,
			bw_est,
			islog,
			0
		)

		out.gamma = beta * A'

		if (dograd) {
			dA = h1_A_matrix(
				delta,
				estimator,
				K,
				cutoff_orig,
				bw_orig,
				cutoff_est,
				bw_est,
				islog,
				1
			)

			out.dgamma_dbeta  = A	
			out.dgamma_ddelta = dA * beta'
		}
		else {
			out.dgamma_dbeta  = J(0, 0, .)
			out.dgamma_ddelta = J(0, 1, .)
		}

		return(out)
	}
				
		// design row transformation for h1, consistent with h1coef_map()
	struct hdesign_out scalar h1design23(
		real scalar delta,
		real colvector zR,
		real scalar cutoff_orig,
		real scalar bw_orig,
		real scalar cutoff_est,
		real scalar bw_est,
		real scalar K,
		real scalar estimator,
		real scalar islog,
		real scalar dograd
	)
	{
		struct hdesign_out scalar out
		real matrix Xbase

		Xbase = pbasis(zR, K)

		out.X = Xbase * h1_A_matrix(
			delta,
			estimator,
			K,
			cutoff_orig,
			bw_orig,
			cutoff_est,
			bw_est,
			islog,
			0
		)

		if (dograd) {
			out.dXddelta = Xbase * h1_A_matrix(
				delta,
				estimator,
				K,
				cutoff_orig,
				bw_orig,
				cutoff_est,
				bw_est,
				islog,
				1
			)
		}
		else {
			out.dXddelta = J(0, 0, .)
		}

		return(out)
	}

		// -----------------------------------------------------------------------------
		// Mass-row helpers for profile estimators 2/3
		// -----------------------------------------------------------------------------

		real rowvector bmodel_row23(
			real scalar delta,
			real scalar cutoff_orig,
			real scalar bw_orig,
			real scalar cutoff_est,
			real scalar bw_est,
			real scalar K,
			real scalar estimator,
			real scalar islog,
			real scalar zbar_est
		)
		{
			real scalar r
			real rowvector R

			if (1 + delta <= 0) {
				_error(3498, "delta must be greater than -1")
			}

			if (estimator == 2) {
				/*
					Chetty restriction. Estimator 2's model is h1 = h0/(1+delta),
					so the mass that left the region above the kink and piled up
					at z* is
						B = int_{zstar}^{zbar} (h0 - h1)
						  = int_{zstar}^{zbar} h0 * (1 - 1/(1+delta))
						  = delta/(1+delta) * int_{zstar}^{zbar} h0(z) dz
						  = delta * int_{zstar}^{zbar} h1(z) dz.
					i.e. delta multiplies the integral of the SHIFTED density h1,
					not the counterfactual h0.  The support endpoints do not move
					(vertical rescaling), so both integrals run to the same zbar.
				*/
				R = (delta / ((1 + delta) * bw_est)) * intbasis(cutoff_est, zbar_est, K)
			}
			else if (estimator == 3) {
				/* Theoretically consistent restriction: bw * B = int_{zstar}^{zstar+r(delta)} h0(z) dz */
				r = response_length(delta, cutoff_orig, bw_orig, bw_est, islog)
				R = intbasis(cutoff_est, cutoff_est + r, K) / bw_est
			}
			else {
				_error(3498, "bmodel_row23 only handles estimators 2 and 3")
			}

			return(R)
		}
		real rowvector d_bmodel_row_ddelta(
		real scalar delta,
		real scalar cutoff_orig,
		real scalar bw_orig,
		real scalar cutoff_est,
		real scalar bw_est,
		real scalar K,
		real scalar estimator,
		real scalar islog,
		real scalar zbar_est,
		real scalar ntheta
	)
	{
		real scalar Kb
		real rowvector out

		out = J(1, ntheta, 0)
		Kb  = K + 1

		if (estimator == 0 | estimator == 1) {
			return(out)
		}

		if (estimator == 2 | estimator == 3) {
			out[1, 1..Kb] = d_bmodel_row23_ddelta(
				delta,
				cutoff_orig,
				bw_orig,
				cutoff_est,
				bw_est,
				K,
				estimator,
				islog,
				zbar_est
			)

			return(out)
		}

		_error(3498, "d_bmodel_row_ddelta only handles estimators 0, 1, 2, and 3")
	}


	real rowvector d_bmodel_row23_ddelta(
		real scalar delta,
		real scalar cutoff_orig,
		real scalar bw_orig,
		real scalar cutoff_est,
		real scalar bw_est,
		real scalar K,
		real scalar estimator,
		real scalar islog,
		real scalar zbar_est
	)
	{
		real scalar r, dr

		if (1 + delta <= 0) {
			_error(3498, "delta must be greater than -1")
		}

		if (estimator == 2) {
			/* d/ddelta [ delta/(1+delta) ] = 1/(1+delta)^2 */
			return(intbasis(cutoff_est, zbar_est, K) / ((1 + delta)^2 * bw_est))
		}

		if (estimator == 3) {
			r  = response_length(delta, cutoff_orig, bw_orig, bw_est, islog)
			dr = d_response_length_ddelta(delta, cutoff_orig, bw_orig, bw_est, islog)

			return(pbasis_row(cutoff_est + r, K) * dr / bw_est)
		}

		_error(3498, "d_bmodel_row23_ddelta only handles estimators 2 and 3")
	}

		
		// -----------------------------------------------------------------------------
		// Unified stacked design/profile objective for estimators 0/1/2/3
		// -----------------------------------------------------------------------------

		real rowvector cf_mass_row(
			real scalar delta,
			real scalar cutoff_orig,
			real scalar bw_orig,
			real scalar cutoff_est,
			real scalar bw_est,
			real scalar K,
			real scalar estimator,
			real scalar islog,
			real scalar zL_excl_orig,
			real scalar zH_excl_orig,
			real scalar zL_excl_est,
			real scalar zH_excl_est,
			real scalar ntheta,
			real scalar nosplit
		)
		{
			real scalar  Kb
			real rowvector R, Rlo, Rhi
			struct hcoef_out scalar h1map

			Kb = K + 1
			R = J(1, ntheta, 0)

			if (estimator == 0) {
				Rlo = intbasis(zL_excl_est, cutoff_est, K) / bw_est
				Rhi = intbasis(cutoff_est, zH_excl_est, K) / bw_est
				R[1, 1..Kb] = Rlo
				R[1, (Kb+1)..(2*Kb)] = Rhi
			}
			else if (estimator == 1) {
				R[1, 1..Kb] = intbasis(zL_excl_est, zH_excl_est, K) / bw_est
			}
			else if (estimator == 2) {
				/*
					Excluded-region counterfactual entering the Chetty
					mass-balance row.
					  splitmass: h0 below z*, h1 = h0/(1+delta) above z*, so
					    (with the bmodel_row term) the restriction is
					    M_E - int_{zL}^{z*} h0 - int_{z*}^{zH} h1
					        = delta/(1+delta) int_{z*}^{zbar} h0.
					  poolmass:  the whole excluded region is counted at h0
					    (no 1/(1+delta) deflation), so the restriction is
					    Chetty's own
					    M_E - int_{zL}^{zH} h0 = delta/(1+delta) int_{z*}^{zbar} h0
					    -- missing mass over the FULL [z*, zbar], and
					    consistent with bunch_transform's poolmass B.
				*/
				Rlo = intbasis(zL_excl_est, cutoff_est, K) / bw_est
				if (nosplit) Rhi = intbasis(cutoff_est, zH_excl_est, K) / bw_est
				else         Rhi = intbasis(cutoff_est, zH_excl_est, K) / ((1 + delta) * bw_est)
				R[1, 1..Kb] = Rlo + Rhi
			}
			else if (estimator == 3) {
				/*
					splitmass: h0 below z*, relocated h1 above z* inside the
					excluded region.  poolmass: whole excluded region at h0,
					so (with the bmodel_row response-interval term) the
					restriction is M_E - int_{zL}^{zH} h0 = int_{z*}^{z*+r} h0,
					matching bunch_transform's poolmass B + eresp() inversion.
				*/
				Rlo = intbasis(zL_excl_est, cutoff_est, K) / bw_est
				if (nosplit) {
					Rhi = intbasis(cutoff_est, zH_excl_est, K) / bw_est
				}
				else {
					h1map = h1coef_map(
						J(1, K+1, 0),
						delta,
						estimator,
						K,
						cutoff_orig,
						bw_orig,
						cutoff_est,
						bw_est,
						islog,
						1
					)
					Rhi = (intbasis(cutoff_est, zH_excl_est, K) * h1map.dgamma_dbeta) / bw_est
				}
				R[1, 1..Kb] = Rlo + Rhi
			}
			else {
				_error(3498, "cf_mass_row only handles estimators 0, 1, 2, and 3")
			}

			return(R)
		}
		
	real rowvector d_cf_mass_row_ddelta(
		real scalar delta,
		real scalar cutoff_orig,
		real scalar bw_orig,
		real scalar cutoff_est,
		real scalar bw_est,
		real scalar K,
		real scalar estimator,
		real scalar islog,
		real scalar zL_excl_orig,
		real scalar zH_excl_orig,
		real scalar zL_excl_est,
		real scalar zH_excl_est,
		real scalar ntheta,
		real scalar nosplit
	)
	{
		real scalar Kb
		real rowvector out, Ihi
		real matrix dA

		out = J(1, ntheta, 0)
		Kb  = K + 1

		if (estimator == 0 | estimator == 1) {
			return(out)
		}

		if (1 + delta <= 0) {
			_error(3498, "delta must be greater than -1")
		}

		/*
			Under poolmass the excluded-region counterfactual in the mass
			row (cf_mass_row) is the plain int_{zL}^{zH} h0 -- no delta -- so
			its delta-derivative is zero for both estimators.  Only the
			splitmass cf row carries delta (through the 1/(1+delta) deflation
			for est 2 and the relocation Jacobian for est 3).
		*/
		if (nosplit) {
			return(out)
		}

		Ihi = intbasis(cutoff_est,zH_excl_est, K)

		if (estimator == 2) {
			out[1, 1..Kb] = -Ihi / ((1 + delta)^2 * bw_est)
		}
		else if (estimator == 3) {
			dA = h1_A_matrix(
				delta,
				estimator,
				K,
				cutoff_orig,
				bw_orig,
				cutoff_est,
				bw_est,
				islog,
				1
			)

			out[1, 1..Kb] = (Ihi * dA) / bw_est
		}
		else {
			_error(3498, "d_cf_mass_row_ddelta only handles estimators 0, 1, 2, and 3")
		}

		return(out)
	}

	real rowvector bmodel_row(
			real scalar delta,
			real scalar cutoff_orig,
			real scalar bw_orig,
			real scalar cutoff_est,
			real scalar bw_est,
			real scalar K,
			real scalar estimator,
			real scalar islog,
			real scalar zbar_est,
			real scalar ntheta
		)
		{
			real scalar r
			real rowvector R

			if (1 + delta <= 0) {
				_error(3498, "delta must be greater than -1")
			}

			R = J(1, ntheta, 0)

			if (estimator == 0 | estimator == 1) {
				return(R)
			}

			if (estimator == 2) {
				/*
					h1 = h0/(1+delta) => missing mass above the kink is
					delta/(1+delta) * int_{zstar}^{zbar} h0  (= delta * int h1).
					See bmodel_row23() for the full derivation.
				*/
				R[1, 1..(K+1)] = (delta / ((1 + delta) * bw_est)) * intbasis(cutoff_est, zbar_est, K)
			}
			else if (estimator == 3) {
				r = response_length(delta, cutoff_orig, bw_orig, bw_est, islog)
				R[1, 1..(K+1)] = intbasis(cutoff_est, cutoff_est + r, K) / bw_est
			}
			else {
				_error(3498, "bmodel_row only handles estimators 0, 1, 2, and 3")
			}

			return(R)
		}



		real colvector make_ystack(
			real colvector y,
			real colvector side,
			real colvector bunch,
			real scalar estimator,
			real scalar Hstar_obs
		)
		{
			if (estimator == 1) {
				return(select(y, bunch :== 0) \ Hstar_obs)
			}

			return(
				select(y, (bunch :== 0) :& (side :== -1)) \
				select(y, (bunch :== 0) :& (side :==  1)) \
				Hstar_obs
			)
		}

		struct design_out scalar make_design(
			real scalar delta,
			real colvector z,
			real colvector side,
			real colvector bunch,
			real scalar cutoff_orig,
			real scalar bw_orig,
			real scalar cutoff_est,
			real scalar bw_est,
			real scalar K,
			real scalar estimator,
			real scalar islog,
			real scalar zL_excl_orig,
			real scalar zH_excl_orig,
			real scalar zL_excl_est,
			real scalar zH_excl_est,
			real scalar zbar_est,
			real scalar dograd,
			real scalar nosplit
		)
		{
			struct design_out scalar out
			struct hdesign_out scalar h1

			real scalar Kb, nL, nR, n0, ntheta
			real colvector zL, zR, z0
			real matrix XL, XR, X0, dX0
			real rowvector Xcf, Xbmod, Xmass, dXcf, dXbmod, dXmass

			Kb = K + 1

			if (estimator == 0) {
				zL = select(z, (bunch :== 0) :& (side :== -1))
				zR = select(z, (bunch :== 0) :& (side :==  1))
				XL = pbasis(zL, K)
				XR = pbasis(zR, K)
				nL = rows(XL)
				nR = rows(XR)
				ntheta = 2*Kb + 1

				Xcf = cf_mass_row(delta, cutoff_orig, bw_orig, cutoff_est, bw_est, K, estimator, islog, zL_excl_orig, zH_excl_orig, zL_excl_est, zH_excl_est, ntheta, nosplit)
				Xmass = Xcf
				Xmass[1, ntheta] = 1

				out.X =
					(XL,              J(nL, Kb, 0), J(nL, 1, 0)) \
					(J(nR, Kb, 0),    XR,           J(nR, 1, 0)) \
					Xmass

				out.dXddelta = J(rows(out.X), 0, .)
			}
			else if (estimator == 1) {
				z0 = select(z, bunch :== 0)
				X0 = pbasis(z0, K)
				n0 = rows(X0)
				ntheta = Kb + 1

				Xcf = cf_mass_row(delta, cutoff_orig, bw_orig, cutoff_est, bw_est, K, estimator, islog, zL_excl_orig, zH_excl_orig, zL_excl_est, zH_excl_est, ntheta, nosplit)
				Xmass = Xcf
				Xmass[1, ntheta] = 1

				out.X =
					(X0, J(n0, 1, 0)) \
					Xmass

				out.dXddelta = J(rows(out.X), 0, .)
			}
			else if (estimator == 2 | estimator == 3) {
				zL = select(z, (bunch :== 0) :& (side :== -1))
				zR = select(z, (bunch :== 0) :& (side :==  1))

				XL = pbasis(zL, K)
				h1 = h1design23(delta, zR, cutoff_orig, bw_orig, cutoff_est, bw_est, K, estimator, islog, dograd)

				ntheta = Kb
				Xcf    = cf_mass_row(delta, cutoff_orig, bw_orig, cutoff_est, bw_est, K, estimator, islog, zL_excl_orig, zH_excl_orig, zL_excl_est, zH_excl_est, ntheta, nosplit)
				Xbmod  = bmodel_row(delta, cutoff_orig, bw_orig, cutoff_est, bw_est, K, estimator, islog, zbar_est, ntheta)
				Xmass  = Xcf + Xbmod

				out.X = XL \ h1.X \ Xmass

				if (dograd) {
					dX0    = J(rows(XL), ntheta, 0)
					dXcf   = d_cf_mass_row_ddelta(
						delta,
						cutoff_orig,
						bw_orig,
						cutoff_est,
						bw_est,
						K,
						estimator,
						islog,
						zL_excl_orig,
						zH_excl_orig,
						zL_excl_est,
						zH_excl_est,
						ntheta,
						nosplit
					)

					dXbmod = d_bmodel_row_ddelta(
						delta,
						cutoff_orig,
						bw_orig,
						cutoff_est,
						bw_est,
						K,
						estimator,
						islog,
						zbar_est,
						ntheta
					)

					dXmass = dXcf + dXbmod

					out.dXddelta = dX0 \ h1.dXddelta \ dXmass
				}
				else {
					out.dXddelta = J(0, 0, .)
				}
			}
			else {
				_error(3498, "make_design only handles estimators 0, 1, 2, and 3")
			}

			return(out)
		}

		real rowvector profTheta(
			real scalar delta,
			real colvector y,
			real colvector z,
			real colvector side,
			real colvector bunch,
			real scalar Hstar_obs,
			real scalar cutoff_orig,
			real scalar bw_orig,
			real scalar cutoff_est,
			real scalar bw_est,
			real scalar K,
			real scalar estimator,
			real scalar islog,
			real scalar zL_excl_orig,
			real scalar zH_excl_orig,
			real scalar zL_excl_est,
			real scalar zH_excl_est,
			real scalar zbar_est,
			real scalar nosplit
		)
		{
			real colvector ystack, theta
			struct design_out scalar D

			ystack = make_ystack(y, side, bunch, estimator, Hstar_obs)

			D = make_design(delta, z, side, bunch, cutoff_orig, bw_orig, cutoff_est, bw_est, K, estimator, islog, zL_excl_orig, zH_excl_orig, zL_excl_est, zH_excl_est, zbar_est, 0, nosplit)

			theta = qrsolve(D.X, ystack)
			return(theta')
		}
		real scalar profQ(
			real scalar delta,
			real colvector y,
			real colvector z,
			real colvector side,
			real colvector bunch,
			real scalar Hstar_obs,
			real scalar cutoff_orig,
			real scalar bw_orig,
			real scalar cutoff_est,
			real scalar bw_est,
			real scalar K,
			real scalar estimator,
			real scalar islog,
			real scalar zL_excl_orig,
			real scalar zH_excl_orig,
			real scalar zL_excl_est,
			real scalar zH_excl_est,
			real scalar zbar_est,
			real scalar positive,
			real scalar nosplit
		)
		{
			real colvector ystack, theta, resid
			struct design_out scalar D

			if (estimator == 0 | estimator == 1) {
				return(0)
			}

			if (positive == 1) {
				if (delta <= 0) return(1e300)
			}
			else {
				if (1 + delta <= 1e-8) return(1e300)
			}

			ystack = make_ystack(y, side, bunch, estimator, Hstar_obs)

			D = make_design(
				delta,
				z,
				side,
				bunch,
				cutoff_orig,
				bw_orig,
				cutoff_est,
				bw_est,
				K,
				estimator,
				islog,
				zL_excl_orig,
				zH_excl_orig,
				zL_excl_est,
				zH_excl_est,
				zbar_est,
				0,
				nosplit
			)

			theta = qrsolve(D.X, ystack)
			resid = ystack - D.X * theta

			return(quadcross(resid, resid))
		}


		// -----------------------------------------------------------------------------
		// Variance estimation
		// -----------------------------------------------------------------------------


	/*
		Pearson overdispersion estimate from the reference-bin fit:
		    phi-hat = (1/(n-p)) * sum_j (y_j - mu_j)^2 / mu_j
		Returns 1 when it cannot be formed.  Used for e(dispersion),
		for scale(x2), and (in variance_robust) to lift the
		bunching-mass row off the pure-Poisson floor.
	*/
	real scalar pearson_phi(
		real colvector y_ref,
		real colvector mu_ref,
		real scalar p
	)
	{
		real scalar n
		real colvector m

		n = rows(y_ref)
		if (n <= p | n < 1) return(1)
		m = mu_ref
		if (missing(m) | missing(y_ref)) return(1)
		m = m + (m :<= 0) :* 1e-8
		return( sum((y_ref - mu_ref):^2 :/ m) / (n - p) )
	}

	/* -------------------------------------------------------------------
	   Reference-region goodness of fit for the fitted counterfactual.
	   y, mu are the observed and fitted bin counts on the bins that enter
	   the fit (h0 below the kink, h1(delta-hat) above); the mass-restriction
	   row is NOT passed in.  p is the number of fitted parameters (K+1 for
	   estimator 1, K+2 for 2/3, 3 for Saez).  Poisson / quasi-Poisson
	   family -- the histogram is multinomial(N, .) and, conditional on N,
	   independent Poisson.  Returns, in order:
	     1 nbins  2 np  3 df(=nbins-np)
	     4 deviance      5 pearson_x2
	     6 loglik        7 deviance_null   8 loglik_null   (null: mu_j = ybar)
	     9 r2_dev (Cameron-Windmeijer, 1 - dev/dev_null)
	    10 aic   11 bic   12 qaic   13 qaicc              (quasi: /phi-hat)
	   All missing if the fit region is degenerate or any mu_j <= 0.
	   ------------------------------------------------------------------- */
	real rowvector pbx_gof(real colvector y, real colvector mu, real scalar p)
	{
		real scalar nb, df, dev, x2, ll, ybar, dev0, ll0, phi, r2d, qaic
		real colvector r, pos

		nb = rows(y)
		if (nb < 2 | nb != rows(mu)) return(J(1, 13, .))
		if (missing(y) | missing(mu) | colmin(mu) <= 0) return(J(1, 13, .))

		r   = y :- mu
		x2  = colsum(r:^2 :/ mu)
		pos = selectindex(y :> 0)                         /* 0*log(0) := 0 */
		if (rows(pos) > 0)
			dev = 2 * (colsum(y[pos] :* log(y[pos] :/ mu[pos])) - colsum(r))
		else
			dev = -2 * colsum(r)
		ll = colsum(y :* log(mu) :- mu :- lngamma(y :+ 1))

		ybar = mean(y)
		if (ybar <= 0) return(J(1, 13, .))
		if (rows(pos) > 0) dev0 = 2 * colsum(y[pos] :* log(y[pos] :/ ybar))
		else               dev0 = 0
		ll0 = colsum(y :* log(ybar) :- ybar :- lngamma(y :+ 1))

		df  = nb - p
		phi = pearson_phi(y, mu, p)                       /* == e(dispersion) */
		if (phi < 1 | phi >= .) phi = 1
		r2d  = (dev0 > 0 ? 1 - dev / dev0 : .)
		qaic = -2*ll/phi + 2*p

		return((nb, p, df, dev, x2, ll, dev0, ll0, r2d,
		        -2*ll + 2*p,
		        -2*ll + p*log(nb),
		        qaic,
		        (nb - p - 1 > 0 ? qaic + 2*p*(p+1)/(nb - p - 1) : .)))
	}


	// -----------------------------------------------------------------------------
	// Residual / wild bootstrap over bins  (vce(bootstrap, residual | wild))
	// -----------------------------------------------------------------------------

	/* mean-0, variance-1 multiplier weights: 1 Rademacher, 2 Mammen, 3 Webb */
	real colvector pbx_wild_weights(real scalar n, real scalar code)
	{
		real colvector u, idx, vals
		real scalar pm, a, bb

		u = runiform(n, 1)
		if (code == 2) {
			pm = (sqrt(5) + 1) / (2 * sqrt(5))
			a  = -(sqrt(5) - 1) / 2
			bb =  (sqrt(5) + 1) / 2
			return( (u :< pm) :* a :+ (u :>= pm) :* bb )
		}
		else if (code == 3) {
			vals = (-sqrt(1.5) \ -1 \ -sqrt(0.5) \ sqrt(0.5) \ 1 \ sqrt(1.5))
			idx  = ceil(u :* 6)
			idx  = idx :+ (idx :< 1) :- (idx :> 6)
			return( vals[idx] )
		}
		return( 2 :* (u :>= 0.5) :- 1 )                     /* Rademacher */
	}

	/* Fit the counterfactual polynomial (separately per side) and cache the
	   fitted counts, the reference residuals and the resampling pool. */
	void pbx_resboot_setup(
		string scalar yv, string scalar zv, string scalar sv, string scalar bv,
		real scalar K
	)
	{
		external real colvector pbx_Y
		external real colvector pbx_MU
		external real colvector pbx_RESID
		external real colvector pbx_POOL
		external real colvector pbx_REF
		external real colvector pbx_EX
		external real scalar    pbx_PHI

		real colvector Y, Z, S, Bn, zc, il, ir
		real matrix X, Xl, Xr

		Y  = st_data(., yv)
		Z  = st_data(., zv)
		S  = st_data(., sv)
		Bn = st_data(., bv)

		zc = Z :- mean(Z)
		if (max(abs(zc)) > 0) zc = zc :/ max(abs(zc))
		X = pbasis(zc, K)

		il = selectindex((Bn :== 0) :& (S :== -1))
		ir = selectindex((Bn :== 0) :& (S :==  1))

		pbx_MU = J(rows(Y), 1, .)

		if (rows(il) > cols(X)) {
			Xl = X[il, .]
			pbx_MU[il] = Xl * qrsolve(Xl, Y[il])
		}
		else {
			pbx_MU[il] = Y[il]
		}

		if (rows(ir) > cols(X)) {
			Xr = X[ir, .]
			pbx_MU[ir] = Xr * qrsolve(Xr, Y[ir])
		}
		else {
			pbx_MU[ir] = Y[ir]
		}

		pbx_Y     = Y
		pbx_REF   = il \ ir
		pbx_EX    = selectindex(Bn :!= 0)
		pbx_RESID = J(rows(Y), 1, 0)
		pbx_RESID[pbx_REF] = Y[pbx_REF] :- pbx_MU[pbx_REF]
		pbx_POOL  = pbx_RESID[pbx_REF]

		/* Pearson overdispersion of the reference bins -- used to scale the
		   perturbation on the (much taller) excluded/bunching bins, so the
		   residual/wild bootstrap of the mass matches variance_robust's
		   phi-hat * Hstar mass-row treatment. */
		pbx_PHI = pearson_phi(Y[pbx_REF], pbx_MU[pbx_REF], cols(X))
		if (pbx_PHI >= . | pbx_PHI < 1) pbx_PHI = 1
	}

	/* -----------------------------------------------------------------------
	   Preferred residual-pool builder: use the fitted counterfactual of the
	   MAIN model, m(theta-hat) -- h0 below the kink, h1 above -- passed in as
	   the full-length Stata matrix r_mu_bins (set by profile_run / saez_run,
	   missing off the reference bins).  The bootstrap then resamples
	   y_j - m(theta-hat) around m(theta-hat) and re-runs the same estimator,
	   so the resampling DGP matches the model being estimated.  Excluded bins
	   are handled exactly as in pbx_resboot_setup (sqrt(phi-hat * y_j), the
	   resampling twin of variance_robust's mass row).  Sets r_rbsetup_ok to 1
	   on success; the caller falls back to pbx_resboot_setup (per-side OLS)
	   otherwise.
	   ----------------------------------------------------------------------- */
	void pbx_resboot_setup_model(string scalar yv, string scalar bunchv,
		real scalar phi_in)
	{
		external real colvector pbx_Y
		external real colvector pbx_MU
		external real colvector pbx_RESID
		external real colvector pbx_POOL
		external real colvector pbx_REF
		external real colvector pbx_EX
		external real scalar    pbx_PHI

		real colvector Y, Bn, MB
		real matrix    MBm
		real scalar    phi

		st_numscalar("r_rbsetup_ok", 0)

		MBm = st_matrix("r_mu_bins")
		if (rows(MBm) < 2 | cols(MBm) < 1) return
		MB = MBm[., 1]

		Y  = st_data(., yv)
		if (rows(MB) != rows(Y)) return
		Bn = st_data(., bunchv)

		pbx_REF = selectindex(MB :< .)
		if (rows(pbx_REF) < 2) return
		if (hasmissing(Y[pbx_REF])) return

		pbx_Y     = Y
		pbx_MU    = MB
		pbx_RESID = J(rows(Y), 1, 0)
		pbx_RESID[pbx_REF] = Y[pbx_REF] :- MB[pbx_REF]
		pbx_POOL  = pbx_RESID[pbx_REF]
		pbx_EX    = selectindex(Bn :!= 0)

		phi = phi_in
		if (phi >= . | phi < 1) phi = 1
		pbx_PHI = phi

		st_numscalar("r_rbsetup_ok", 1)
	}

	/* One resampled count vector, written straight back to `yv'. */
	void pbx_resboot_draw(string scalar yv, real scalar wildflag, real scalar wildcode)
	{
		external real colvector pbx_Y
		external real colvector pbx_MU
		external real colvector pbx_RESID
		external real colvector pbx_POOL
		external real colvector pbx_REF
		external real colvector pbx_EX
		external real scalar    pbx_PHI

		real colvector col, ev, w, idx, wex
		real scalar n, npool, nex

		n     = rows(pbx_Y)
		npool = rows(pbx_POOL)
		nex   = rows(pbx_EX)

		col = pbx_Y
		ev  = J(n, 1, 0)

		/* Reference bins: residual = iid pool draw; wild = own residual x weight. */
		if (wildflag == 0) {
			idx = ceil(runiform(npool, 1) :* npool)
			idx = idx :+ (idx :< 1)
			ev[pbx_REF] = pbx_POOL[idx]
		}
		else {
			w = pbx_wild_weights(n, wildcode)
			ev[pbx_REF] = pbx_RESID[pbx_REF] :* w[pbx_REF]
		}

		/* Excluded / bunching bins: perturb with overdispersed-Poisson noise
		   scaled to each bin's own height, sqrt(phi-hat * y_j) * weight --
		   matching variance_robust's phi-hat * Hstar mass row (the reference
		   residual pool is on the wrong scale for the much taller mass bins). */
		if (nex > 0) {
			if (wildflag == 0) wex = rnormal(nex, 1, 0, 1)
			else               wex = pbx_wild_weights(nex, wildcode)
			ev[pbx_EX] = sqrt(pbx_PHI :* abs(pbx_Y[pbx_EX])) :* wex
		}

		col[pbx_REF] = pbx_MU[pbx_REF] :+ ev[pbx_REF]
		if (nex > 0) col[pbx_EX] = pbx_Y[pbx_EX] :+ ev[pbx_EX]
		col = col :* (col :> 0)

		st_store(., yv, col)
	}


	real matrix variance_multinomial(
		real matrix G_stack,
		real colvector y_stack,
		real scalar addcons
	)
	{
		real scalar N
		real matrix G, Vm, bread, V

		G = G_stack

		if (addcons == 1) {
			G = G, J(rows(G), 1, 1)
		}

		if (rows(G) != rows(y_stack)) {
			return(J(cols(G), cols(G), .))
		}

		if (missing(G) | missing(y_stack)) {
			return(J(cols(G), cols(G), .))
		}

		N = sum(y_stack)

		if (N <= 0 | N >= .) {
			return(J(cols(G), cols(G), .))
		}

		Vm = diag(y_stack) - (y_stack * y_stack') / N
		Vm = (Vm + Vm') / 2

		bread = pinv(quadcross(G, G))

		V = bread * G' * Vm * G * bread
		V = (V + V') / 2

		return(V)
	}


	/*
		Residual / Eicker-White delta-method variance of the stacked
		count regression.

		Identical plumbing to variance_multinomial -- same bread
		pinv(G'G), same tiny multinomial off-diagonals -yi*yj/N -- but the
		DIAGONAL of the meat matrix for the polynomial-fit rows is replaced
		by the squared fit residual (y_j - mu_j)^2, with an HC finite-sample
		correction:

		    hctype 0 : (y_j - mu_j)^2                       (HC0)
		    hctype 1 : * n/(n-p)                            (HC1)
		    hctype 2 : / (1 - h_jj)                         (HC2)
		    hctype 3 : / (1 - h_jj)^2                       (HC3)

		with h_jj the leverage g_j (G'G)^-1 g_j'.  n counts only the
		robust (bin) rows and p = cols(G).

		The trailing `nmassrows' row(s) -- the observed bunching-mass
		identity Hstar = int h0 (+ int h1) + B -- keep the multinomial
		variance Hstar*(1 - Hstar/N): that row is an accounting identity
		with no lack of fit, and its residual is ~0 by construction, so a
		squared-residual meat there would spuriously zero out the (large)
		sampling noise in the observed bunching mass.

		Robust to: arbitrary bin-level heteroskedasticity; polynomial
		misspecification treated as exchangeable noise (round-number
		heaping, secondary bumps, a neighbouring kink); the marginal,
		per-bin part of overdispersion from repeated individuals or random
		year effects.

		NOT robust to: cross-bin correlation of any kind (off-diagonals
		stay at the multinomial value) -- individuals migrating between
		adjacent bins across years, common shocks that move a band of bins
		together, the covariance half of a panel design effect.  Cluster on
		the microdata outside polbunch for that.
	*/
	real matrix variance_robust(
		real matrix G_stack,
		real colvector y_stack,
		real colvector mu_stack,
		real scalar hctype,
		real scalar nmassrows,
		real scalar masscorr
	)
	{
		real scalar N, n, p, j, dfadj, hv, phi
		real matrix G, Vm, bread, V, PG
		real colvector e, dm, dr, hjj

		G = G_stack

		if (rows(G) != rows(y_stack) | rows(G) != rows(mu_stack)) {
			return(J(cols(G), cols(G), .))
		}

		if (missing(G) | missing(y_stack) | missing(mu_stack)) {
			return(J(cols(G), cols(G), .))
		}

		N = sum(y_stack)

		if (N <= 0 | N >= .) {
			return(J(cols(G), cols(G), .))
		}

		n = rows(y_stack) - nmassrows
		p = cols(G)

		if (n < 1) {
			return(J(cols(G), cols(G), .))
		}

		Vm = diag(y_stack) - (y_stack * y_stack') / N
		Vm = (Vm + Vm') / 2

		e  = y_stack - mu_stack
		dm = diagonal(Vm)
		dr = e :^ 2

		if (hctype == 1) {
			dfadj = (n > p ? n / (n - p) : 1)
			dr = dr :* dfadj
		}
		else if (hctype == 2 | hctype == 3) {
			bread = pinv(quadcross(G, G))
			PG = G * bread
			hjj = rowsum(PG :* G)
			for (j = 1; j <= rows(hjj); j++) {
				hv = hjj[j]
				if (hv >= 0.9999) hv = 0.9999
				if (hv < 0) hv = 0
				if (hctype == 2) dr[j] = dr[j] / (1 - hv)
				else             dr[j] = dr[j] / (1 - hv)^2
			}
		}

		/*
			Trailing accounting row(s) -- the observed bunching mass.
			Its fitted residual is ~0 by construction, so it cannot get a
			squared-residual meat.  Default (masscorr==0): keep the pure
			multinomial variance.  masscorr==1: lift it by the Pearson
			overdispersion of the reference bins, phi-hat * Hstar(1-Hstar/N)
			-- borrowing the reference-region dispersion into the excluded
			window, which is what Chetty's residual bootstrap does
			implicitly and what closes the estimator-1 gap.
		*/
		phi = 1
		if (masscorr == 1) {
			phi = pearson_phi(y_stack[1::n], mu_stack[1::n], p)
			if (phi >= .) phi = 1
			if (phi < 1)  phi = 1
		}
		for (j = n + 1; j <= rows(y_stack); j++) {
			dr[j] = phi * dm[j]
		}

		if (missing(dr)) {
			return(J(cols(G), cols(G), .))
		}

		Vm = Vm - diag(dm) + diag(dr)
		Vm = (Vm + Vm') / 2

		bread = pinv(quadcross(G, G))

		V = bread * G' * Vm * G * bread
		V = (V + V') / 2

		/*
			Swapping in the squared-residual diagonal above while keeping
			the (unchanged) multinomial off-diagonals is not guaranteed to
			leave V positive semi-definite -- when the residual-implied
			diagonal is small relative to the retained off-diagonal
			covariances, V can pick up a genuinely negative eigenvalue
			(confirmed empirically: as large in magnitude as other real
			variances, not floating-point noise). A non-PSD theta-level V
			propagates through ANY subsequent delta-method transform
			(bunch_transform's G V G' is a congruence, which preserves
			indefiniteness), and Stata's own -ereturn post- silently
			substitutes an all-zero V -- with only a terse, easily-missed
			warning -- when it judges a posted V "nonsymmetric or highly
			singular".  Repair by projecting V onto the PSD cone: clip
			negative eigenvalues to 0 and reconstruct.  Flagged via
			r_robust_psdclip so the .ado layer can expose e(vce_psdclip)
			and warn once on the main fit.
		*/
		{
			real matrix Ue
			real rowvector lambda
			real scalar negtol

			symeigensystem(V, Ue, lambda)
			negtol = -1e-8 * max((max(abs(lambda)), 1))
			if (min(lambda) < negtol) {
				lambda = lambda :* (lambda :> 0)
				V = Ue * diag(lambda) * Ue'
				V = (V + V') / 2
				st_numscalar("r_robust_psdclip", 1)
			}
			else {
				st_numscalar("r_robust_psdclip", 0)
			}
		}

		return(V)
	}


	// -----------------------------------------------------------------------------
	// Cluster-robust variance from binned data  (vce(cluster clustvar))
	// -----------------------------------------------------------------------------

	/* Stash the individual (cluster id, bin id) pairs before polbunch
	   collapses the data to a histogram. */
	void pbx_cluster_stash(string scalar idv, string scalar binv)
	{
		external real colvector pbx_CLID
		external real colvector pbx_CLBIN

		pbx_CLID  = st_data(., idv)
		pbx_CLBIN = st_data(., binv)
	}

	/* Build M_raw = sum_i c_i c_i', the raw J x J bin co-visitation matrix
	   over all J histogram bins (bin order), where c_i is cluster i's
	   vector of visit counts across the bins.  binposname is a B x 2 Stata
	   matrix mapping each bin id (col 1) to its position 1..J (col 2). */
	void pbx_covis_build(string scalar binposname, real scalar nbin)
	{
		external real colvector pbx_CLID
		external real colvector pbx_CLBIN

		real matrix bp, A, M
		real colvector lut, id, pos, cv, nz
		real scalar mn, r, i, s, e, a, b, n

		bp = st_matrix(binposname)
		mn = min(bp[., 1])
		lut = J(max(bp[., 1]) - mn + 1, 1, .)
		for (r = 1; r <= rows(bp); r++) {
			lut[bp[r, 1] - mn + 1] = bp[r, 2]
		}

		id  = pbx_CLID
		pos = lut[pbx_CLBIN :- (mn - 1)]
		pos = editmissing(pos, 0)

		A = (id, pos)
		_sort(A, 1)
		id  = A[., 1]
		pos = A[., 2]

		M = J(nbin, nbin, 0)
		n = 0
		i = 1
		while (i <= rows(id)) {
			s = i
			while (i <= rows(id)) {
				if (id[i] != id[s]) break
				i++
			}
			e = i - 1

			cv = J(nbin, 1, 0)
			for (r = s; r <= e; r++) {
				if (pos[r] >= 1 & pos[r] <= nbin) cv[pos[r]] = cv[pos[r]] + 1
			}
			nz = selectindex(cv)
			for (a = 1; a <= rows(nz); a++) {
				for (b = 1; b <= rows(nz); b++) {
					M[nz[a], nz[b]] = M[nz[a], nz[b]] + cv[nz[a]] * cv[nz[b]]
				}
			}
			n++
		}

		st_matrix("r_pbx_Mraw", M)
		st_numscalar("r_pbx_nclust", n)
	}

	/* Fold the raw J x J co-visitation matrix onto the stacked estimating
	   system.  srowname is a J x 1 Stata matrix sending each raw bin to its
	   stacked row (reference bins keep their row; every excluded bin -> the
	   trailing mass row).  Result r_pbx_Mstack is Jstack x Jstack. */
	void pbx_covis_collapse(string scalar Mrawname, string scalar srowname, real scalar Jstack)
	{
		real matrix Mraw, M
		real colvector sr
		real scalar nb, i, j, a, b

		Mraw = st_matrix(Mrawname)
		sr   = st_matrix(srowname)
		if (cols(sr) > 1) sr = sr[., 1]
		nb   = rows(Mraw)

		M = J(Jstack, Jstack, 0)
		for (i = 1; i <= nb; i++) {
			a = sr[i]
			if (a < 1 | a > Jstack) continue
			for (j = 1; j <= nb; j++) {
				b = sr[j]
				if (b < 1 | b > Jstack) continue
				M[a, b] = M[a, b] + Mraw[i, j]
			}
		}
		M = (M + M') / 2
		st_matrix("r_pbx_Mstack", M)
	}

	/* Validate a user-supplied co-visitation matrix against the histogram.
	   r_pbx_covchk: 0 ok, 1 not symmetric, 2 diag(M) < y somewhere,
	   3 wrong dimension. */
	void pbx_covis_check(string scalar Mname, string scalar yname)
	{
		real matrix M
		real colvector y, d
		real scalar tol

		M = st_matrix(Mname)
		y = st_data(., yname)

		st_numscalar("r_pbx_covchk", 0)

		if (rows(M) != cols(M) | rows(M) != rows(y)) {
			st_numscalar("r_pbx_covchk", 3)
			return
		}
		tol = 1e-6 * (1 + max(abs(vec(M))))
		if (max(abs(vec(M - M'))) > tol) {
			st_numscalar("r_pbx_covchk", 1)
			return
		}
		d = diagonal(M) - y
		if (min(d) < -1e-6 * (1 + max(y))) {
			st_numscalar("r_pbx_covchk", 2)
			return
		}
	}

	/* Cluster-robust sandwich: meat = n/(n-1) * (M - y y' / n). */
	real matrix variance_cluster(
		real matrix G_stack,
		real colvector y_stack,
		real matrix M,
		real scalar nclust
	)
	{
		real matrix G, meat, bread, V

		G = G_stack

		if (rows(G) != rows(y_stack) | rows(M) != rows(y_stack) | cols(M) != rows(y_stack)) {
			return(J(cols(G), cols(G), .))
		}
		if (missing(G) | missing(y_stack) | missing(M) | nclust < 2 | nclust >= .) {
			return(J(cols(G), cols(G), .))
		}

		meat  = (nclust / (nclust - 1)) :* (M - (y_stack * y_stack') / nclust)
		meat  = (meat + meat') / 2
		bread = pinv(quadcross(G, G))
		V     = bread * G' * meat * G * bread
		V     = (V + V') / 2

		return(V)
	}



		// -----------------------------------------------------------------------------
		// Bunching-response inversion and transformed output
		// -----------------------------------------------------------------------------
		real scalar eresp(
			real scalar B,
			real scalar tau,
			real rowvector cf,
			real scalar bw,
			real scalar xscale
		)
		{
			real scalar target, lo, hi, mid
			real scalar Flo, Fhi, Fmid
			real scalar iter, maxiter
			real rowvector cfpoly, intpoly

			/*
				cf is [b1, ..., bK, b0].
				poly* uses [b0, b1, ..., bK].
			*/
			if (cols(cf) == 1) {
				cfpoly = cf
			}
			else {
				cfpoly = cf[cols(cf)], cf[1..cols(cf)-1]
			}

			target = B * bw

			if (target <= 0 | target >= .) {
				return(.)
			}

			/*
				Constant case.
			*/
			if (cols(cfpoly) == 1) {
				if (cfpoly[1] <= 0 | cfpoly[1] >= .) {
					return(.)
				}
				return((target / cfpoly[1]) * xscale)
			}

			intpoly = polyinteg(cfpoly, 1)

			/*
				F(r_est) = integral_tau^{tau+r_est} h0(u)du - target.
				eresp returns r in original units, so final multiply by xscale.
			*/
			lo  = 0
			Flo = -target

			hi  = 1
			Fhi = polyeval(intpoly, tau + hi) - polyeval(intpoly, tau) - target

			/*
				Expand upper bracket until crossing, but cap to avoid infinite search.
				In normalized case, hi is in normalized/bin units.
				In non-normalized case, hi is in original z units.
			*/
			while (Fhi <= 0 & hi < 1e6) {
				hi  = 2 * hi
				Fhi = polyeval(intpoly, tau + hi) - polyeval(intpoly, tau) - target
			}

			if (Fhi <= 0 | Fhi >= .) {
				return(.)
			}

			maxiter = 100

			for (iter = 1; iter <= maxiter; iter++) {
				mid  = (lo + hi) / 2
				Fmid = polyeval(intpoly, tau + mid) - polyeval(intpoly, tau) - target

				if (Fmid >= 0) {
					hi = mid
				}
				else {
					lo = mid
				}
			}

			return(hi * xscale)
		}



		void bunch_transform(
			real rowvector theta,
			real scalar estimator,
			real scalar K,
			real scalar cutoff_orig,
			real scalar cutoff_est,
			real scalar bw_orig,
			real scalar bw_est,
			real scalar xscale,
			real scalar islog,
			real scalar constant,
			real scalar hastax,
			real scalar t0,
			real scalar t1,
			real scalar dograd,
			real scalar zbar_est,
			real scalar Hstar_obs,
			real scalar nosplit,
			real scalar zL_excl_orig,
			real scalar zH_excl_orig,
			real scalar zL_excl_est,
			real scalar zH_excl_est,
			real scalar weakid
		)
		{
			struct hcoef_out scalar h1

			real rowvector xcut, dm, dB, RB, dRB
			real rowvector beta, gamma
			real rowvector dr, dF_dbeta, dMR, dshift
			real matrix G
			real rowvector b

			real scalar Kb, nout, i, hasresp, e3delta
			real rowvector Rlo, Rhi, Ihi
			real scalar m, EM, B, delta
			real scalar r, u, h_u, shift, MR, elast, A

			real rowvector ibeta, igamma
			real scalar idelta, iBraw
			real rowvector obeta, ogamma
			real scalar oB, oEM, odelta, oshift, oMR, oe

			// PARSE RAW PARAMS
			Kb = K + 1

			ibeta = 1..Kb
			beta  = theta[ibeta]

			if (estimator == 0) {
				igamma = (Kb+1)..(2*Kb)
				iBraw  = 2*Kb + 1
				gamma  = theta[igamma]
				B      = theta[iBraw]
			}
			else if (estimator == 1) {
				iBraw = Kb + 1
				gamma = beta
				B     = theta[iBraw]
			}
			else {
				idelta = Kb + 1
				delta  = theta[idelta]
			}

			// ALLOCATE OUTPUT
			hasresp = 1

			/*
				A structural-delta column is reported for estimator 2
				always, and for estimator 3 when the reported response is
				NOT the pure structural delta -- under constant (constant-
				density inversion) or under poolmass (nosplit==1, response
				backed out from the pooled reduced-form B).
			*/
			e3delta = (estimator==3 & (constant | nosplit))

			nout =
				2*Kb +                                      // h0,h1
				2 +                                         // B, EM
				(estimator==2 | e3delta) +                  // delta
				2 +                                         // shift, MR
				hastax

			b = J(1, nout, .)
			if (dograd) G = J(nout, cols(theta), 0)

			// output indices
			obeta  = ibeta
			ogamma = (Kb+1)..(2*Kb)
			oB     = 2*Kb + 1
			oEM    = 2*Kb + 2

			if (estimator==2 | e3delta) {
				i = 1
				odelta = 2*Kb + 3
			}
			else i = 0

			oshift = 2*Kb + 3 + i
			oMR    = 2*Kb + 4 + i
			if (hastax) oe = 2*Kb + 5 + i

			// h0
			b[1,obeta] = beta
			if (dograd) G[obeta,ibeta] = I(Kb)

			// h1
			if (estimator == 0) {
				b[1,ogamma] = gamma
				if (dograd) G[ogamma,igamma] = I(Kb)
			}
			else if (estimator == 1) {
				b[1,ogamma] = beta
				if (dograd) G[ogamma,ibeta] = I(Kb)
			}
			else {
			   h1 = h1coef_map(
					beta,
					delta,
					estimator,
					K,
					cutoff_orig,
					bw_orig,
					cutoff_est,
					bw_est,
					islog,
					dograd
				)

				b[1,ogamma] = h1.gamma

				if (dograd) {
					G[ogamma,ibeta]  = h1.dgamma_dbeta
					G[ogamma,idelta] = h1.dgamma_ddelta
				}
			}

			// B
			//   poolmass  (nosplit==1): B = M - int_{zL}^{zH} h0hat
			//   splitmass (nosplit==0): B = M - int_{zL}^{z*} h0hat
			//                                 - int_{z*}^{zH} h1hat
			// Estimator 1 has h0==h1, so the two coincide -- keep the free
			// profile B and ignore nosplit.
			if (estimator == 1) {
				b[1,oB] = B
				if (dograd) G[oB,iBraw] = 1
			}
			else if (estimator == 0) {
				if (nosplit) {
					// poolmass: extrapolate the fitted left polynomial h0hat
					// across the whole excluded region
					RB = intbasis(zL_excl_est, zH_excl_est, K) / bw_est
					B  = Hstar_obs - RB * beta'
					if (dograd) G[oB, ibeta] = -RB
				}
				else {
					// splitmass: the stacked-fit B row already equals
					//   M - int_{zL}^{z*} h0hat - int_{z*}^{zH} h1hat
					// (h1hat is the separately fitted right polynomial)
					if (dograd) G[oB, iBraw] = 1
				}
				b[1,oB] = B
			}
			else if (estimator == 2) {
				if (nosplit) {
					// poolmass (default): reduced-form B from the pooled
					// h0hat integral (Chetty-style)
					RB = intbasis(zL_excl_est, zH_excl_est, K) / bw_est
					B  = Hstar_obs - RB * beta'
					if (dograd) {
						G[oB, ibeta]  = -RB
						G[oB, idelta] = 0
					}
				}
				else {
					// splitmass: h1hat = h0hat/(1+delta) on the right
					Rlo = intbasis(zL_excl_est, cutoff_est, K) / bw_est
					Ihi = intbasis(cutoff_est, zH_excl_est, K)
					Rhi = Ihi / ((1 + delta) * bw_est)
					B   = Hstar_obs - (Rlo + Rhi) * beta'
					if (dograd) {
						G[oB, ibeta]  = -(Rlo + Rhi)
						G[oB, idelta] = (Ihi / ((1 + delta)^2 * bw_est)) * beta'
					}
				}
				b[1,oB] = B
			}
			else if (estimator == 3) {
				if (nosplit) {
					// poolmass: reduced-form B from the pooled h0hat integral
					RB = intbasis(zL_excl_est, zH_excl_est, K) / bw_est
					B  = Hstar_obs - RB * beta'
					if (dograd) {
						G[oB, ibeta]  = -RB
						G[oB, idelta] = 0
					}
				}
				else {
					// splitmass (default): theoretically consistent
					// model-implied B,  bw*B = int_{z*}^{z*+r(delta)} h0hat
					RB = bmodel_row23(delta, cutoff_orig, bw_orig, cutoff_est, bw_est, K, estimator, islog, zbar_est)
					B  = RB * beta'
					if (dograd) {
						G[oB, ibeta] = RB
						dRB = d_bmodel_row23_ddelta(delta, cutoff_orig, bw_orig, cutoff_est, bw_est, K, estimator, islog, zbar_est)
						G[oB, idelta] = dRB * beta'
					}
				}
				b[1,oB] = B
			}

			if (dograd) dB = G[oB,.]

			// Excess mass
			xcut = pbasis_row(cutoff_est, K)
			m = beta * xcut'

			EM = B / m
			b[1,oEM] = EM

			if (dograd) {
				dm = J(1, cols(theta), 0)
				dm[ibeta] = xcut
				G[oEM,.] = dB/m - (B/(m^2))*dm
			}

			// structural delta column (see e3delta above)
			if (estimator == 2 | e3delta) {
				b[1,odelta] = delta
				if (dograd) G[odelta,idelta] = 1
			}

			/*
				Estimator 2/3 profile delta was weakly identified (boundary of
				the deltamax() search, or a multi-modal concentrated SSR).
				Report the counterfactual, the structural delta column, B and
				excess mass, but drop the response: shift / marginal_response /
				elasticity are backed out of delta through a steep nonlinear
				map and are not credible here.  Same truncation as the
				no-real-root path below.
			*/
			if (weakid) {
				hasresp = 0
				b = b[1, 1..(oshift-1)]
				if (dograd) G = G[1..(oshift-1), .]

				st_numscalar("b_bunchcalc_hasresp", hasresp)
				st_matrix("b_bunchcalc", b)
				if (dograd) {
					st_matrix("G_bunchcalc", G)
				}
				else {
					st_matrix("G_bunchcalc", J(0,0,.))
				}
				return
			}

			// Shift, marginal response, elasticity
			if (hastax) A = ln(1-t0) - ln(1-t1)

			if (constant) {
				MR = B*bw_orig/m

				if (islog == 0) {
					shift = MR/cutoff_orig
				}
				else {
					shift = exp(MR) - 1
				}

				b[1,oshift] = shift
				b[1,oMR]    = MR

				if (dograd) {
					dMR = bw_orig * (dB/m - (B/(m^2))*dm)

					if (islog == 0) {
						dshift = dMR/cutoff_orig
					}
					else {
						dshift = exp(MR)*dMR
					}

					G[oshift,.] = dshift
					G[oMR,.]    = dMR
				}
			}
			else if (estimator == 3 & !nosplit) {
				// estimator 3, splitmass: shift is the structural delta
				// (model-consistent bunching mass; exact inversion)
				if (islog == 0) {
					shift = delta
					MR    = delta*cutoff_orig
				}
				else {
					shift = delta
					MR    = ln(1+delta)
				}

				b[1,oshift] = shift
				b[1,oMR]    = MR

				if (dograd) {
					G[oshift,idelta] = 1

					if (islog == 0) {
						G[oMR,idelta] = cutoff_orig
					}
					else {
						G[oMR,idelta] = 1/(1+delta)
					}
				}
			}
			else {
				// estimators 0/1/2, and estimator 3 under poolmass:
				// back the response out of B by solving the counterfactual-
				// density integral equation with eresp()
				r = eresp(B, cutoff_est, beta, bw_est, xscale)

				if (r >= .) {
					/*
						If no admissible root exists for the missing-mass
						equation, shift/marginal_response/elasticity are
						not identified for this draw -- drop them from the
						posted vector (as before). An earlier attempt to
						post them as ordinary Stata missing values instead
						ran into "ereturn post"/"ereturn repost" both
						unconditionally refusing any b vector containing
						missing entries ("matrix has missing values",
						r(504), with no override option in this Stata
						version) -- so the width genuinely has to vary
						with hasresp. The .ado layer now exposes
						e(hasresp) instead, so callers can check that
						BEFORE touching _b[bunching:elasticity] rather
						than needing capture at all.
					*/
					hasresp = 0
					b = b[1, 1..(oshift-1)]
					if (dograd) G = G[1..(oshift-1), .]

					st_numscalar("b_bunchcalc_hasresp", hasresp)
					st_matrix("b_bunchcalc", b)
					if (dograd) {
						st_matrix("G_bunchcalc", G)
					}
					else {
						st_matrix("G_bunchcalc", J(0,0,.))
					}
					return
				}

				u = cutoff_est + r/xscale
				h_u = beta * pbasis_row(u,K)'

				b[1,oMR] = r

				if (islog == 0) {
					shift = r/cutoff_orig
				}
				else {
					shift = exp(r) - 1
				}

				b[1,oshift] = shift

				if (dograd) {
					dF_dbeta = intbasis(cutoff_est, u, K)

					dr = (xscale*bw_est/h_u)*dB
					dr[ibeta] = dr[ibeta] - (xscale/h_u)*dF_dbeta

					G[oMR,.] = dr

					if (islog == 0) {
						G[oshift,.] = dr/cutoff_orig
					}
					else {
						G[oshift,.] = exp(r)*dr
					}
				}
			}

			if (hastax) {
				elast = ln(1+shift)/A
				b[1,oe] = elast

				if (dograd) {
					G[oe,.] = G[oshift,.] / ((1+shift)*A)
				}
			}

			st_numscalar("b_bunchcalc_hasresp", hasresp)
			st_matrix("b_bunchcalc", b)

			if (dograd) {
				st_matrix("G_bunchcalc", G)
			}
			else {
				st_matrix("G_bunchcalc", J(0,0,.))
			}
		}

		// -----------------------------------------------------------------------------
		// Top-level unified profile estimator
		// -----------------------------------------------------------------------------

		/*
			Multi-start solve of the estimator 2/3 profile over the scalar
			delta.  The concentrated SSR profQ(delta) is minimised by a coarse
			candidate grid, then a golden-section refinement around EVERY
			sampled interior local minimum -- the same shape as
			polbunch_mdt_mata's delta search.

			This replaces a single Newton start (optimize()), which lands in
			different basins depending on initdelta / polynomial order: when
			the level-shift restriction is misspecified the profile SSR is
			multi-modal and near-singular in places, so the starting value
			decides the answer.

			delta is confined to [dlo, dhi].  profQ() already returns 1e300
			for an infeasible delta, so the bounded search cannot throw.

			Written back by reference:
			  bestQ   SSR at the returned optimum
			  nbasin  number of interior local minima found on the grid
			  gapQ    (2nd-best basin SSR)/(best basin SSR), or . if nbasin<2
		*/
		real scalar prof_delta_solve(
			real colvector y,
			real colvector z,
			real colvector side,
			real colvector bunch,
			real scalar Hstar_obs,
			real scalar cutoff_orig,
			real scalar bw_orig,
			real scalar cutoff_est,
			real scalar bw_est,
			real scalar K,
			real scalar estimator,
			real scalar islog,
			real scalar zL_excl_orig,
			real scalar zH_excl_orig,
			real scalar zL_excl_est,
			real scalar zH_excl_est,
			real scalar zbar_est,
			real scalar positive,
			real scalar nosplit,
			real scalar initdelta,
			real scalar dlo,
			real scalar dhi,
			real scalar bestQ,
			real scalar nbasin,
			real scalar gapQ
		)
		{
			real colvector grid, cand, gd, gq, rds, rqs
			real scalar lo, hi, j, d, Q, ng, gr
			real scalar a, b, x1, x2, f1, f2, iter
			real scalar rd, rQ, bd, bQ, b2Q, btol, dgap

			bestQ  = .
			nbasin = 0
			gapQ   = .

			lo = dlo
			hi = dhi
			if (hi <= lo) return(.)

			if (positive == 1) {
				grid = (1e-8 \ 1e-6 \ 1e-4 \ 1e-3 \ 0.005 \ 0.01 \ 0.025 \
				        0.05 \ 0.075 \ 0.10 \ 0.15 \ 0.20 \ 0.30 \ 0.50 \
				        0.75 \ 1 \ 1.5 \ 2 \ 3 \ 5)
			}
			else {
				grid = (-0.99 \ -0.95 \ -0.90 \ -0.80 \ -0.70 \ -0.60 \ -0.50 \
				        -0.40 \ -0.30 \ -0.20 \ -0.15 \ -0.10 \ -0.075 \ -0.05 \
				        -0.025 \ -0.01 \ 0 \ 0.01 \ 0.025 \ 0.05 \ 0.075 \ 0.10 \
				        0.15 \ 0.20 \ 0.30 \ 0.50 \ 0.75 \ 1 \ 1.5 \ 2 \ 3 \ 5)
			}

			if (initdelta > lo & initdelta < hi) {
				cand = sort(grid \ initdelta, 1)
			}
			else {
				cand = sort(grid, 1)
			}

			// evaluate the feasible, de-duplicated candidates once
			gd = J(0, 1, .)
			gq = J(0, 1, .)
			for (j = 1; j <= rows(cand); j++) {
				d = cand[j]
				if (d <= lo | d >= hi) continue
				if (rows(gd) >= 1) {
					if (d == gd[rows(gd)]) continue
				}
				Q = profQ(d, y, z, side, bunch, Hstar_obs, cutoff_orig, bw_orig, cutoff_est, bw_est, K, estimator, islog, zL_excl_orig, zH_excl_orig, zL_excl_est, zH_excl_est, zbar_est, positive, nosplit)
				gd = gd \ d
				gq = gq \ Q
			}

			ng = rows(gd)
			if (ng == 0) return(.)

			gr  = (sqrt(5) - 1) / 2
			rds = J(0, 1, .)
			rqs = J(0, 1, .)

			// golden-section refine around every local minimum of the grid
			for (j = 1; j <= ng; j++) {

				if (j > 1) {
					if (gq[j] > gq[j-1]) continue
				}
				if (j < ng) {
					if (gq[j] > gq[j+1]) continue
				}

				// bracket = neighbouring grid points, widened to the feasible
				// bounds at the ends so a near-boundary optimum is reachable
				if (j > 1)  a = gd[j-1]
				else        a = lo
				if (j < ng) b = gd[j+1]
				else        b = hi
				if (a < lo) a = lo
				if (b > hi) b = hi

				if (b <= a) {
					rd = gd[j]
					rQ = gq[j]
				}
				else {
					x1 = b - gr * (b - a)
					x2 = a + gr * (b - a)
					f1 = profQ(x1, y, z, side, bunch, Hstar_obs, cutoff_orig, bw_orig, cutoff_est, bw_est, K, estimator, islog, zL_excl_orig, zH_excl_orig, zL_excl_est, zH_excl_est, zbar_est, positive, nosplit)
					f2 = profQ(x2, y, z, side, bunch, Hstar_obs, cutoff_orig, bw_orig, cutoff_est, bw_est, K, estimator, islog, zL_excl_orig, zH_excl_orig, zL_excl_est, zH_excl_est, zbar_est, positive, nosplit)

					for (iter = 1; iter <= 100; iter++) {
						if (abs(b - a) < 1e-10 * max((1, abs(x1), abs(x2)))) break
						if (f1 > f2) {
							a  = x1
							x1 = x2
							f1 = f2
							x2 = a + gr * (b - a)
							f2 = profQ(x2, y, z, side, bunch, Hstar_obs, cutoff_orig, bw_orig, cutoff_est, bw_est, K, estimator, islog, zL_excl_orig, zH_excl_orig, zL_excl_est, zH_excl_est, zbar_est, positive, nosplit)
						}
						else {
							b  = x2
							x2 = x1
							f2 = f1
							x1 = b - gr * (b - a)
							f1 = profQ(x1, y, z, side, bunch, Hstar_obs, cutoff_orig, bw_orig, cutoff_est, bw_est, K, estimator, islog, zL_excl_orig, zH_excl_orig, zL_excl_est, zH_excl_est, zbar_est, positive, nosplit)
						}
					}

					if (f1 <= f2) {
						rd = x1
						rQ = f1
					}
					else {
						rd = x2
						rQ = f2
					}

					// never let the refinement beat the sampled point itself
					if (gq[j] < rQ) {
						rd = gd[j]
						rQ = gq[j]
					}
				}

				rds = rds \ rd
				rqs = rqs \ rQ
			}

			// global optimum among the refined minima
			bQ = 1e300
			bd = .
			for (j = 1; j <= rows(rqs); j++) {
				if (rqs[j] < bQ) {
					bQ = rqs[j]
					bd = rds[j]
				}
			}

			// no interior local minimum found: take the best grid point outright
			if (bd >= .) {
				for (j = 1; j <= ng; j++) {
					if (gq[j] < bQ) {
						bQ = gq[j]
						bd = gd[j]
					}
				}
			}

			bestQ = bQ

			/*
				Competing-basin check for the weak-ID guard.  A runner-up
				counts only if it is
				  - a genuine INTERIOR minimum (not the search sliding into the
				    lb / deltamax bound -- a monotone tail truncated by the cap
				    is not a second solution), and
				  - at a materially different delta from the winner.
				nbasin then counts the winner plus every such competitor, and
				gapQ is (best competitor SSR)/(winner SSR).
			*/
			btol = 1e-6
			dgap = 0.05
			if (0.25 * abs(bd) > dgap) dgap = 0.25 * abs(bd)

			nbasin = 1
			b2Q    = 1e300
			for (j = 1; j <= rows(rds); j++) {
				if (rds[j] == bd) continue
				if (rds[j] <= lo + btol | rds[j] >= hi - btol) continue
				if (abs(rds[j] - bd) <= dgap) continue
				nbasin = nbasin + 1
				if (rqs[j] < b2Q) b2Q = rqs[j]
			}
			if (nbasin >= 2 & bQ > 0 & b2Q < 1e299) gapQ = b2Q / bQ

			return(bd)
		}

		void profile_run(
			string scalar yvar,
			string scalar zvar,
			string scalar sidevar,
			string scalar bunchvar,
			real scalar cutoff_orig,
			real scalar bw_orig,
			real scalar cutoff_est,
			real scalar bw_est,
			real scalar K,
			real scalar estimator,
			real scalar islog,
			real scalar zL_excl_orig,
			real scalar zH_excl_orig,
			real scalar zL_excl_est,
			real scalar zH_excl_est,
			real scalar zbar_est,
			real scalar dovar,
			real scalar initdelta,
			real scalar positive,
			real scalar hctype,
			real scalar masscorr,
			real scalar deltamax,
			string scalar covisname,
			real scalar nclust,
			real scalar nosplit

		)
		{
			real colvector y, z, side, bunch
			real scalar Kb, Hstar_obs, delta_hat, lb, nref_
			real scalar dhi, bestQ, nbasin, gapQ, weakid
			real rowvector theta_hat, beta_hat, b
			real matrix Vout

			struct design_out scalar D
			real colvector ystack, mu
			real colvector gof_grs, gof_gib, gof_gia
			real scalar    gof_gnb
			real matrix Gv


			Kb = K + 1

			y     = st_data(., yvar)
			z     = st_data(., zvar)
			side  = st_data(., sidevar)
			bunch = st_data(., bunchvar)

			Hstar_obs = sum(select(y, bunch :> 0))

			if (initdelta <= 0 | initdelta >= .) {
				initdelta = 0.05
			}

			if (estimator == 0 | estimator == 1) {
				delta_hat = 0
				theta_hat = profTheta(
					delta_hat,
					y,
					z,
					side,
					bunch,
					Hstar_obs,
					cutoff_orig,
					bw_orig,
					cutoff_est,
					bw_est,
					K,
					estimator,
					islog,
					zL_excl_orig,
					zH_excl_orig,
					zL_excl_est,
					zH_excl_est,
					zbar_est,
					nosplit
				)
				b = theta_hat
			}
			else if (estimator == 2 | estimator == 3) {

				if (positive == 1) lb = 1e-8
				else               lb = -1 + 1e-8
				if (initdelta <= lb | initdelta >= .) initdelta = 0.05

				/*
					Economic upper bound on the structural delta: deltamax
					(default 1 -- a 100% earnings response at the kink).
				*/
				dhi = deltamax
				if (dhi >= . | dhi <= 0)  dhi = 1
				if (dhi <= lb)            dhi = lb + 1e-6

				bestQ  = .
				nbasin = .
				gapQ   = .
				delta_hat = prof_delta_solve(
					y, z, side, bunch, Hstar_obs,
					cutoff_orig, bw_orig, cutoff_est, bw_est,
					K, estimator, islog,
					zL_excl_orig, zH_excl_orig, zL_excl_est, zH_excl_est,
					zbar_est, positive, nosplit, initdelta, lb, dhi,
					bestQ, nbasin, gapQ)

				if (delta_hat >= .) delta_hat = initdelta

				/*
					Weak identification of delta:
					  - the search hit a bound of [lb, dhi]  (runaway response
					    or a boundary optimum), or
					  - the concentrated SSR has >=2 competing local minima
					    (near-tied basins, ratio < 1.10).
					When flagged, the .ado withholds shift / marginal_response
					/ elasticity (e(hasresp)=0) and reports e(delta_weakid)=1.
				*/
				weakid = 0
				if (delta_hat <= lb + 1e-6)         weakid = 1
				if (delta_hat >= dhi - 1e-6)        weakid = 1
				if (nbasin >= 2 & gapQ < 1.10)      weakid = 1

				st_numscalar("r_weakid_profile", weakid)
				st_numscalar("r_nbasin_profile", nbasin)

				beta_hat = profTheta(
					delta_hat,
					y,
					z,
					side,
					bunch,
					Hstar_obs,
					cutoff_orig,
					bw_orig,
					cutoff_est,
					bw_est,
					K,
					estimator,
					islog,
					zL_excl_orig,
					zH_excl_orig,
					zL_excl_est,
					zH_excl_est,
					zbar_est,
					nosplit
				)

				theta_hat = beta_hat, delta_hat
				b = theta_hat
			}
			else {
				_error(3498, "profile_run only handles estimators 0, 1, 2, and 3")
			}

			if (dovar == 1) {
				D = make_design(delta_hat, z, side, bunch,
					cutoff_orig, bw_orig, cutoff_est, bw_est, K, estimator, islog,
					zL_excl_orig, zH_excl_orig, zL_excl_est, zH_excl_est, zbar_est, 1, nosplit)

				ystack = make_ystack(y, side, bunch, estimator, Hstar_obs)

				if (estimator == 0 | estimator == 1) {
					Gv = D.X
					mu = D.X * theta_hat'
				}
				else {
					beta_hat = theta_hat[1, 1..Kb]
					Gv = D.X, (D.dXddelta * beta_hat')
					mu = D.X * beta_hat'
				}

				if (covisname != "") {
					Vout = variance_cluster(Gv, ystack, st_matrix(covisname), nclust)
				}
				else if (hctype >= 0) {
					Vout = variance_robust(Gv, ystack, mu, hctype, 1, masscorr)
				}
				else {
					Vout = variance_multinomial(Gv, ystack, 0)
				}

				/* Pearson overdispersion phi-hat from the reference rows,
				   returned for e(dispersion) and scale(x2). */
				nref_ = rows(ystack) - 1
				if (nref_ >= 1) {
					st_numscalar("r_phi_profile", pearson_phi(ystack[1::nref_], mu[1::nref_], cols(Gv)))
					st_matrix("r_gof", pbx_gof(ystack[1::nref_], mu[1::nref_], cols(Gv)))
					st_numscalar("r_gof_rmse", sqrt(mean((ystack[1::nref_] :- mu[1::nref_]):^2)))
					st_numscalar("r_gof_massresid",
						(ystack[nref_ + 1] > 0
						 ? (ystack[nref_ + 1] - mu[nref_ + 1]) / ystack[nref_ + 1]
						 : .))

					/*
						Reference-region fit split at the cutoff.  Below z*
						nobody bunches and h1 = h0, so A1/A2 are vacuous:
						overdispersion / deviance there is a near-pure test of
						A3 (the polynomial's fit to h0).  Above z* it also picks
						up an A1/A2 shape error in h1.  df per side uses the h0
						order K+1.  Stacked reference order is [below ; above]
						for estimators 0/2/3, dataset order for estimator 1.
					*/
					gof_grs = select(side, bunch :== 0)
					if (estimator == 1) {
						gof_gib = selectindex(gof_grs :== -1)
						gof_gia = selectindex(gof_grs :==  1)
					}
					else {
						gof_gnb = colsum((bunch :== 0) :& (side :== -1))
						gof_gib = J(0, 1, .)
						if (gof_gnb >= 1)     gof_gib = (1::gof_gnb)
						gof_gia = J(0, 1, .)
						if (gof_gnb < nref_)  gof_gia = ((gof_gnb + 1)::nref_)
					}
					st_matrix("r_gof_below",
						(rows(gof_gib) >= 2 ? pbx_gof(ystack[gof_gib], mu[gof_gib], K + 1) : J(1, 13, .)))
					st_matrix("r_gof_above",
						(rows(gof_gia) >= 2 ? pbx_gof(ystack[gof_gia], mu[gof_gia], K + 1) : J(1, 13, .)))
				}
			}
			else {
				Vout = J(cols(b), cols(b), .)
			}

			st_matrix("r_b_profile", b)
			st_matrix("r_V_profile", Vout)

			if (dovar == 1) {
				st_matrix("r_G_stack", Gv)
				st_matrix("r_mu_stack", mu)
				st_matrix("r_ystack", ystack)

				/*
					Full-length fitted-mean vector for the residual bootstrap:
					m(theta-hat) on the reference bins -- h0(z_j) below the
					kink, h1(z_j; delta-hat) above -- and missing on the
					excluded bins (which carry no per-bin counterfactual, only
					the aggregate mass restriction).  vce(bootstrap, residual|
					wild) resamples y_j - m(theta-hat) around THIS, so the
					bootstrap DGP is the same model that is estimated.
				*/
				{
					real colvector rb_idx, rb_mu
					real scalar    rb_nref

					rb_nref = rows(mu) - 1
					if (estimator == 1) {
						rb_idx = selectindex(bunch :== 0)
					}
					else {
						rb_idx = selectindex((bunch :== 0) :& (side :== -1)) \
						         selectindex((bunch :== 0) :& (side :==  1))
					}
					rb_mu = J(rows(y), 1, .)
					if (rows(rb_idx) == rb_nref & rb_nref >= 1) {
						rb_mu[rb_idx] = mu[1::rb_nref]
					}
					/* st_matrix caps at ~11000 rows; on a finer histogram the
					   residual bootstrap falls back to the per-side fit. */
					if (rows(rb_mu) <= 10000) st_matrix("r_mu_bins", rb_mu)
				}
			}
		}

	real scalar delta_from_mass_e3(
		real rowvector beta0,
		real scalar B,
		real scalar cutoff_orig,
		real scalar bw_orig,
		real scalar cutoff_est,
		real scalar bw_est,
		real scalar K,
		real scalar islog
	)
	{
		real rowvector beta, grid, R
		real scalar d0, d1, f0, f1, mid, fmid
		real scalar i, iter

		beta = beta0
		if (rows(beta) > 1) beta = beta'

		if (cols(beta) != K + 1) {
			_error(3200, "delta_from_mass_e3(): beta has wrong length")
		}

		if (B >= .) return(.)

		grid = (
			-0.999, -0.99, -0.98, -0.95, -0.90, -0.80, -0.70, -0.60,
			-0.50, -0.40, -0.30, -0.20, -0.15, -0.10, -0.075, -0.05,
			-0.025, -0.01, -0.005, 0, 0.005, 0.01, 0.025, 0.05,
			0.075, 0.10, 0.15, 0.20, 0.30, 0.50, 0.75, 1, 1.5, 2,
			3, 5, 10
		)

		d0 = grid[1]
		R  = bmodel_row23(
			d0,
			cutoff_orig,
			bw_orig,
			cutoff_est,
			bw_est,
			K,
			3,
			islog,
			0
		)
		f0 = R * beta' - B

		if (f0 == 0) return(d0)

		for (i = 2; i <= cols(grid); i++) {
			d1 = grid[i]

			R = bmodel_row23(
				d1,
				cutoff_orig,
				bw_orig,
				cutoff_est,
				bw_est,
				K,
				3,
				islog,
				0
			)
			f1 = R * beta' - B

			if (f1 == 0) return(d1)

			if (f0 < . & f1 < . & f0*f1 < 0) {
				for (iter = 1; iter <= 100; iter++) {
					mid = (d0 + d1) / 2

					R = bmodel_row23(
						mid,
						cutoff_orig,
						bw_orig,
						cutoff_est,
						bw_est,
						K,
						3,
						islog,
						0
					)
					fmid = R * beta' - B

					if (fmid == 0) return(mid)

					if (f0*fmid <= 0) {
						d1 = mid
						f1 = fmid
					}
					else {
						d0 = mid
						f0 = fmid
					}
				}

				return((d0 + d1) / 2)
			}

			d0 = d1
			f0 = f1
		}

		return(.)
	}

	real matrix polbunch_mdt_qG(
		real rowvector theta,
		real matrix V,
		real scalar delta,
		real scalar estimator,
		real scalar cutoff_orig,
		real scalar bw_orig,
		real scalar cutoff_est,
		real scalar bw_est,
		real scalar K,
		real scalar islog,
		real scalar zbar_est,
		real scalar positive
	)
	{
		real scalar Kb, ntheta, B
		real rowvector beta, gamma, R
		real colvector q
		real matrix Gq
		struct hcoef_out scalar hmap

		Kb = K + 1
		ntheta = 2*Kb + 1

		if (cols(theta) != ntheta) return(J(0,0,.))
		if (rows(V) != ntheta | cols(V) != ntheta) return(J(0,0,.))
		if (missing(theta) | missing(V)) return(J(0,0,.))

		if (estimator == 1) {
			delta = 0
		}
		else {
			if (delta >= .) return(J(0,0,.))

			if (positive == 1) {
				if (delta <= 1e-8) return(J(0,0,.))
			}
			else {
				if (1 + delta <= 1e-8) return(J(0,0,.))
			}
		}

		beta  = theta[1, 1..Kb]
		gamma = theta[1, (Kb+1)..(2*Kb)]
		B     = theta[1, 2*Kb + 1]

		if (estimator == 1) {
			/*
				q(delta) = gamma - beta
				No nuisance delta.
			*/
			q = (gamma - beta)'

			Gq = J(Kb, ntheta, 0)
			Gq[., 1..Kb]          = -I(Kb)
			Gq[., (Kb+1)..(2*Kb)] =  I(Kb)

			return((q, Gq))
		}

		/*
			For estimators 2 and 3, stack:
				q_shape = gamma - gamma_model(beta, delta)
				q_mass  = B - B_model(beta, delta)

			Important:
				bmodel_row23() is already correct:
				  estimator 2 uses zbar_est;
				  estimator 3 ignores zbar_est and uses the response interval.
		*/
		hmap = h1coef_map(
			beta,
			delta,
			estimator,
			K,
			cutoff_orig,
			bw_orig,
			cutoff_est,
			bw_est,
			islog,
			1
		)

		R = bmodel_row23(
			delta,
			cutoff_orig,
			bw_orig,
			cutoff_est,
			bw_est,
			K,
			estimator,
			islog,
			zbar_est
		)
	
		q = J(Kb + 1, 1, .)
		q[1..Kb, 1] = (gamma - hmap.gamma)'
		q[Kb+1, 1] = B - beta * R'

		Gq = J(Kb + 1, ntheta, 0)

		/*
			Conditional Jacobian with respect to unrestricted theta,
			holding delta fixed. The minimization over delta accounts for
			the one fitted nuisance parameter through df = rows(q) - 1.
		*/
		Gq[1..Kb, 1..Kb]          = -hmap.dgamma_dbeta
		Gq[1..Kb, (Kb+1)..(2*Kb)] =  I(Kb)

		Gq[Kb+1, 1..Kb]     = -R
		Gq[Kb+1, 2*Kb + 1]  =  1

		return((q, Gq))
	}

	real scalar polbunch_mdt_crit(
		real scalar delta,
		real rowvector theta,
		real matrix V,
		real scalar estimator,
		real scalar cutoff_orig,
		real scalar bw_orig,
		real scalar cutoff_est,
		real scalar bw_est,
		real scalar K,
		real scalar islog,
		real scalar zbar_est,
		real scalar positive
	)
	{
		real matrix qG, Gq, Vq, Vqi
		real colvector q
		real scalar m, W

		qG = polbunch_mdt_qG(
			theta,
			V,
			delta,
			estimator,
			cutoff_orig,
			bw_orig,
			cutoff_est,
			bw_est,
			K,
			islog,
			zbar_est,
			positive
		)

		if (rows(qG) == 0) return(1e300)

		m  = rows(qG)
		q  = qG[., 1]
		Gq = qG[., 2..cols(qG)]

		if (missing(q) | missing(Gq)) return(1e300)

		Vq = Gq * V * Gq'
		if (missing(Vq)) return(1e300)

		Vqi = pinv(Vq)
		if (missing(Vqi)) return(1e300)

		W = (q' * Vqi * q)[1,1]

		if (W < 0 | W >= .) return(1e300)

		return(W)
	}

	void polbunch_mdt_mata(
		string scalar bname,
		string scalar Vname,
		real scalar estimator,
		real scalar cutoff_orig,
		real scalar bw_orig,
		real scalar cutoff_est,
		real scalar bw_est,
		real scalar K,
		real scalar islog,
		real scalar zbar_est,
		real scalar positive,
		real scalar initdelta
	)
	{
		real scalar Kb, ntheta, df, pval
		real scalar delta_hat, W_hat
		real scalar lo, span
		real scalar j, d, W, bestd, bestW
		real scalar a, b, x1, x2, f1, f2, gr, iter
		real rowvector theta, grid, candidates
		real matrix V

		Kb = K + 1
		ntheta = 2*Kb + 1

		st_numscalar("r(pb_md)", .)
		st_numscalar("r(pb_md_p)", .)
		st_numscalar("r(pb_md_df)", .)
		st_numscalar("r(pb_md_delta)", .)
		st_numscalar("r(pb_md_failcode)", 0)

		if (!(estimator == 1 | estimator == 2 | estimator == 3)) {
			st_numscalar("r(pb_md_failcode)", 101)
			return
		}

		theta = st_matrix(bname)
		V     = st_matrix(Vname)

		if (cols(theta) != ntheta | rows(V) != ntheta | cols(V) != ntheta) {
			st_numscalar("r(pb_md_failcode)", 102)
			return
		}

		if (missing(theta) | missing(V)) {
			st_numscalar("r(pb_md_failcode)", 103)
			return
		}


		/*
			Estimator 1 has no nuisance delta.
			Test gamma = beta directly.
		*/
		if (estimator == 1) {
			W_hat = polbunch_mdt_crit(
				0,
				theta,
				V,
				estimator,
				cutoff_orig,
				bw_orig,
				cutoff_est,
				bw_est,
				K,
				islog,
				zbar_est,
				positive
			)

			if (W_hat >= 1e299 | W_hat >= .) {
				st_numscalar("r(pb_md_failcode)", 201)
				return
			}

			df = Kb
			pval = chi2tail(df, W_hat)

			st_numscalar("r(pb_md)", W_hat)
			st_numscalar("r(pb_md_p)", pval)
			st_numscalar("r(pb_md_df)", df)
			st_numscalar("r(pb_md_delta)", .)
			st_numscalar("r(pb_md_failcode)", 0)
			return
		}

		/*
			Estimators 2 and 3: minimize over scalar delta.

			The search is deliberately simple and robust:
			  1. Build a legal candidate grid including initdelta.
			  2. Pick the best grid point.
			  3. Golden-section refine between its neighboring grid points.
		*/

		if (positive == 1) {
			lo = 1e-8
			if (initdelta <= lo | initdelta >= .) initdelta = 0.05

			grid = (
				1e-8, 1e-6, 1e-4, 1e-3, 0.005, 0.01, 0.025,
				0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.50,
				0.75, 1, 1.5, 2
			)
		}
		else {
			lo = -1 + 1e-8
			if (initdelta <= lo | initdelta >= .) initdelta = 0.05

			grid = (
				-0.999999, -0.999, -0.99, -0.98, -0.95, -0.90,
				-0.80, -0.70, -0.60, -0.50, -0.40, -0.30,
				-0.20, -0.15, -0.10, -0.075, -0.05, -0.025,
				-0.01, -0.005, 0, 0.005, 0.01, 0.025, 0.05,
				0.075, 0.10, 0.15, 0.20, 0.30, 0.50, 0.75,
				1, 1.5, 2
			)
		}

		/*
			Add initdelta explicitly. Sort manually by evaluating all candidates;
			no need to physically sort for the coarse step.
		*/
		candidates = grid, initdelta

		bestd = .
		bestW = 1e300

		for (j = 1; j <= cols(candidates); j++) {
			d = candidates[j]

			if (d <= lo) continue

			W = polbunch_mdt_crit(
				d,
				theta,
				V,
				estimator,
				cutoff_orig,
				bw_orig,
				cutoff_est,
				bw_est,
				K,
				islog,
				zbar_est,
				positive
			)

			if (W < bestW) {
				bestW = W
				bestd = d
			}
		}

		if (bestd >= . | bestW >= 1e299) {
			st_numscalar("r(pb_md_failcode)", 202)
			return
		}

		/*
			Local bracket around bestd.

			This is deliberately conservative. The criterion is one-dimensional;
			even if it is not globally convex, the coarse grid chooses a sensible
			basin and the local refinement improves the minimum inside that basin.
		*/
		if (bestd > 0) {
			a = max((lo, bestd / 2))
			b = bestd * 2
		}
		else {
			span = max((0.05, abs(bestd) / 2))
			a = max((lo, bestd - span))
			b = bestd + span
		}

		if (b <= a) {
			a = max((lo, bestd - 0.05))
			b = bestd + 0.05
		}

		/*
			Golden-section minimization over [a,b].
		*/
		gr = (sqrt(5) - 1) / 2

		x1 = b - gr * (b - a)
		x2 = a + gr * (b - a)

		f1 = polbunch_mdt_crit(
			x1,
			theta,
			V,
			estimator,
			cutoff_orig,
			bw_orig,
			cutoff_est,
			bw_est,
			K,
			islog,
			zbar_est,
			positive
		)

		f2 = polbunch_mdt_crit(
			x2,
			theta,
			V,
			estimator,
			cutoff_orig,
			bw_orig,
			cutoff_est,
			bw_est,
			K,
			islog,
			zbar_est,
			positive
		)

		for (iter = 1; iter <= 100; iter++) {
			if (abs(b - a) < 1e-10 * max((1, abs(x1), abs(x2)))) break

			if (f1 > f2) {
				a  = x1
				x1 = x2
				f1 = f2
				x2 = a + gr * (b - a)

				f2 = polbunch_mdt_crit(
					x2,
					theta,
					V,
					estimator,
					cutoff_orig,
					bw_orig,
					cutoff_est,
					bw_est,
					K,
					islog,
					zbar_est,
					positive
				)
			}
			else {
				b  = x2
				x2 = x1
				f2 = f1
				x1 = b - gr * (b - a)

				f1 = polbunch_mdt_crit(
					x1,
					theta,
					V,
					estimator,
					cutoff_orig,
					bw_orig,
					cutoff_est,
					bw_est,
					K,
					islog,
					zbar_est,
					positive
				)
			}
		}

		if (f1 <= f2) {
			delta_hat = x1
			W_hat = f1
		}
		else {
			delta_hat = x2
			W_hat = f2
		}

		/*
			Do not let the local refinement make things worse than the coarse grid.
		*/
		if (bestW < W_hat) {
			delta_hat = bestd
			W_hat = bestW
		}

		if (W_hat >= 1e299 | W_hat >= . | delta_hat >= .) {
			st_numscalar("r(pb_md_failcode)", 203)
			return
		}

		/*
			rows(q) = Kb + 1 for estimators 2/3.
			We minimized over one nuisance scalar delta.
			df = Kb.
		*/
		df = Kb
		pval = chi2tail(df, W_hat)

		st_numscalar("r(pb_md)", W_hat)
		st_numscalar("r(pb_md_p)", pval)
		st_numscalar("r(pb_md_df)", df)
		st_numscalar("r(pb_md_delta)", delta_hat)
		st_numscalar("r(pb_md_failcode)", 0)
	}

	/*
		polbunch_shausman_mata -- scalar Hausman test on the ELASTICITY.

		Contrasts the restricted estimator's elasticity e_R (efficient under
		its own H0) with e_U, the elasticity implied by the unrestricted
		two-sided fit (estimator 0).  The minimum-distance test is an omnibus
		chi2_{K+1} on the whole cross-kink coefficient vector; this is a
		focused chi2_1 on the single parameter of interest.

		The two elasticities come from the SAME bin counts and covary, so the
		denominator is the variance of the CONTRAST, built from the joint
		influence functions:

			e_j - e0  ~=  psi_j' (y - mu),   psi_j = A_j' grad_j'

		A_j = diag(1/s_j) pinv(G_j / s_j) is the stacked-fit "bread" and
		grad_j is the delta-method gradient of the elasticity w.r.t. the raw
		stacked coefficients -- the LAST row of the bunch_transform Jacobian
		e(G).  Then

			Var(e_R - e_U) = (psi_R - psi_U)' Vm (psi_R - psi_U)

		with Vm the multinomial (or cluster / HC) meat.  No efficiency
		assumption on e_R is needed -- this is the generalised (Wooldridge)
		Hausman form, valid for Chetty as well as the model-consistent
		estimator.  Both e_R and e_U use A1's marginal-buncher inversion to
		map excess mass to an elasticity, so the test targets A1's density
		restriction given A3, NOT the shared inversion premise.

		When e_U is weakly identified (wild extrapolation of h0 into the
		excluded region) the denominator inflates and the test loses power
		rather than over-rejecting; the returned seU lets the caller show it.
	*/
	void polbunch_shausman_mata(
		string scalar GbcUname,
		string scalar GUname,
		real scalar eU,
		string scalar GbcRname,
		string scalar GRname,
		real scalar eR,
		string scalar yname,
		string scalar muUname,
		real scalar hctype,
		string scalar covUname,
		real scalar nclust
	)
	{
		real scalar N, stat, varD, seU, nbinU, pU, jj, adjU
		real colvector y, muU, evec, dmU, drU, psiU, psiR, psiD, hxU
		real rowvector gradU, gradR, sU, sR
		real matrix GbcU, GbcR, GU, GR, Vm, GUs, GRs, AU, AR, Mclu, PGU

		st_numscalar("r(pb_sh_chi2)", .)
		st_numscalar("r(pb_sh_p)", .)
		st_numscalar("r(pb_sh_df)", .)
		st_numscalar("r(pb_sh_seU)", .)
		st_numscalar("r(pb_sh_failcode)", 0)

		if (missing(eU) | missing(eR)) {
			st_numscalar("r(pb_sh_failcode)", 110)
			return
		}

		GbcU = st_matrix(GbcUname)
		GbcR = st_matrix(GbcRname)
		GU   = st_matrix(GUname)
		GR   = st_matrix(GRname)
		y    = st_matrix(yname)
		if (rows(y) == 1) y = y'

		if (missing(GbcU) | missing(GbcR) | missing(GU) | missing(GR) | missing(y)) {
			st_numscalar("r(pb_sh_failcode)", 101)
			return
		}
		if (rows(GU) != rows(y) | rows(GR) != rows(y)) {
			st_numscalar("r(pb_sh_failcode)", 102)
			return
		}

		/* elasticity gradient = last row of each transform Jacobian */
		gradU = GbcU[rows(GbcU), .]
		gradR = GbcR[rows(GbcR), .]
		if (cols(gradU) != cols(GU) | cols(gradR) != cols(GR)) {
			st_numscalar("r(pb_sh_failcode)", 103)
			return
		}
		if (missing(gradU) | missing(gradR)) {
			st_numscalar("r(pb_sh_failcode)", 111)
			return
		}

		N = sum(y)
		if (N <= 0 | N >= .) {
			st_numscalar("r(pb_sh_failcode)", 104)
			return
		}

		/* ---- meat: multinomial, then cluster / HC as requested ---- */
		Vm = diag(y) - (y * y') / N
		Vm = (Vm + Vm') / 2

		if (covUname != "" & covUname != "." & nclust >= 2) {
			Mclu = st_matrix(covUname)
			if (rows(Mclu) == rows(y) & cols(Mclu) == rows(y) & !missing(Mclu)) {
				Vm = (nclust / (nclust - 1)) :* (Mclu - (y * y') / nclust)
				Vm = (Vm + Vm') / 2
			}
		}

		if (hctype >= 0 & muUname != "" & muUname != ".") {
			muU = st_matrix(muUname)
			if (rows(muU) == 1) muU = muU'
			if (rows(muU) == rows(y) & !missing(muU)) {
				evec  = y - muU
				dmU   = diagonal(Vm)
				drU   = evec :^ 2
				nbinU = rows(y) - 1
				pU    = cols(GU)
				if (hctype == 1) {
					adjU = (nbinU > pU ? nbinU / (nbinU - pU) : 1)
					drU  = drU :* adjU
				}
				else if (hctype == 2 | hctype == 3) {
					PGU = GU * pinv(quadcross(GU, GU))
					hxU = rowsum(PGU :* GU)
					for (jj = 1; jj <= nbinU; jj++) {
						if (hxU[jj] < 1) {
							if (hctype == 2) drU[jj] = drU[jj] / (1 - hxU[jj])
							else             drU[jj] = drU[jj] / (1 - hxU[jj])^2
						}
					}
				}
				for (jj = nbinU + 1; jj <= rows(y); jj++) drU[jj] = dmU[jj]
				if (!missing(drU)) {
					Vm = Vm - diag(dmU) + diag(drU)
					Vm = (Vm + Vm') / 2
				}
			}
		}

		/* ---- influence functions of the two elasticities ---- */
		sU = sqrt(colsum(GU:^2))
		sR = sqrt(colsum(GR:^2))
		if (any(sU :<= 0) | any(sU :>= .) | any(sR :<= 0) | any(sR :>= .)) {
			st_numscalar("r(pb_sh_failcode)", 105)
			return
		}
		GUs = GU :/ sU
		GRs = GR :/ sR
		AU = diag(1 :/ sU) * pinv(GUs)
		AR = diag(1 :/ sR) * pinv(GRs)

		psiU = AU' * gradU'
		psiR = AR' * gradR'
		psiD = psiR - psiU

		varD = (psiD' * Vm * psiD)[1, 1]
		seU  = sqrt((psiU' * Vm * psiU)[1, 1])

		if (varD <= 0 | varD >= .) {
			st_numscalar("r(pb_sh_failcode)", 106)
			return
		}

		stat = (eR - eU)^2 / varD

		st_numscalar("r(pb_sh_chi2)", stat)
		st_numscalar("r(pb_sh_p)", chi2tail(1, stat))
		st_numscalar("r(pb_sh_df)", 1)
		st_numscalar("r(pb_sh_seU)", seU)
		st_numscalar("r(pb_sh_failcode)", 0)
	}



		void polbunch_wald_from_unrestricted(
			string scalar bname,
			string scalar Vname,
			real scalar estimator,
			real scalar cutoff_orig,
			real scalar bw_orig,
			real scalar cutoff_est,
			real scalar bw_est,
			real scalar K,
			real scalar islog,
			real scalar zbar_est
		)
		{
			real scalar Kb, B, delta, Rbeta, Fdelta
			real scalar W, pval, df

			real rowvector theta, beta, gamma
			real rowvector R, Rbmod
			real rowvector ddelta_dbeta, ddelta_dtheta

			real matrix V, Gq, Vq
			real colvector q

			struct hcoef_out scalar hmap

			Kb = K + 1

			theta = st_matrix(bname)
			V     = st_matrix(Vname)

			beta  = theta[1, 1..Kb]
			gamma = theta[1, (Kb+1)..(2*Kb)]
			B     = theta[1, 2*Kb + 1]

			/*
				Initialize outputs as missing.
			*/
			st_numscalar("r(pb_wald)", .)
			st_numscalar("r(pb_p)", .)
			st_numscalar("r(pb_df)", .)
			st_numscalar("r(pb_delta_U)", .)

			if (estimator == 1) {
				/*
					Estimator 1 restrictions:
						gamma = beta

					q = gamma - beta
				*/
				q = (gamma - beta)'

				Gq = J(Kb, 2*Kb + 1, 0)
				Gq[., 1..Kb]           = -I(Kb)
				Gq[., (Kb+1)..(2*Kb)]  =  I(Kb)

				delta = .
			}
			else if (estimator == 2) {
				/*
					Estimator 2:
						gamma = beta / (1 + delta)
						B = delta/(1+delta) * int_{cutoff}^{zbar} h0(z) dz / bw
						  = delta/(1+delta) * Rbeta
					(the missing-mass term is delta * int h1, not delta * int h0,
					because h1 = h0/(1+delta) -- see bmodel_row23()).

					Invert the mass equation for delta_U:
						B (1+delta) = delta Rbeta
						=> delta_U = B / (Rbeta - B)
				*/

				R = intbasis(cutoff_est, zbar_est, K) / bw_est
				Rbeta = R * beta'

				if (Rbeta <= 0 | Rbeta >= . | B >= . | Rbeta - B <= 0) return

				delta = B / (Rbeta - B)

				if (1 + delta <= 0 | delta >= .) return

				q = (gamma - beta :/ (1 + delta))'

				/*
					delta = B / (Rbeta - B),  Rbeta = R * beta'
					ddelta/dbeta_j = -B * R_j / (Rbeta - B)^2
					ddelta/dB      =  Rbeta   / (Rbeta - B)^2
				*/
				ddelta_dbeta = -B * R / ((Rbeta - B)^2)

				ddelta_dtheta =
					ddelta_dbeta,
					J(1, Kb, 0),
					Rbeta / ((Rbeta - B)^2)

				/*
					q_j = gamma_j - beta_j/(1+delta)

					dq_j/dbeta =
						-e_j/(1+delta)
						+ beta_j/(1+delta)^2 * ddelta/dbeta

					dq_j/dgamma = e_j

					dq_j/dB =
						beta_j/(1+delta)^2 * ddelta/dB
				*/
				Gq = J(Kb, 2*Kb + 1, 0)

				Gq[., 1..Kb] =
					-I(Kb)/(1 + delta) +
					(beta' * ddelta_dtheta[1, 1..Kb]) / ((1 + delta)^2)

				Gq[., (Kb+1)..(2*Kb)] = I(Kb)

				Gq[., 2*Kb + 1] =
					beta' * ddelta_dtheta[1, 2*Kb + 1] / ((1 + delta)^2)
			}
			else if (estimator == 3) {
				/*
					Estimator 3:
						B = Bmodel(beta, delta)
						gamma = gamma_map(beta, delta)

					Use the unrestricted mass equation to define delta_U implicitly:

						F(beta, B, delta) = beta * R_B(delta)' - B = 0

					where:
						R_B(delta) = bmodel_row23(delta, ..., estimator=3)

					Then:
						F_beta  = R_B(delta)
						F_B     = -1
						F_delta = dR_B(delta)/ddelta * beta'

					so:
						ddelta/dbeta = -F_beta/F_delta
						ddelta/dB    =  1/F_delta
				*/

				delta = delta_from_mass_e3(
					beta,
					B,
					cutoff_orig,
					bw_orig,
					cutoff_est,
					bw_est,
					K,
					islog
				)

				if (delta >= .) {
					st_numscalar("r(pb_failcode)", 301)
					return
				}

				if (1 + delta <= 0) {
					st_numscalar("r(pb_failcode)", 302)
					return
				}

				hmap = h1coef_map(
					beta,
					delta,
					3,
					K,
					cutoff_orig,
					bw_orig,
					cutoff_est,
					bw_est,
					islog,
					1
				)

				q = (gamma - hmap.gamma)'

				Rbmod = bmodel_row23(
					delta,
					cutoff_orig,
					bw_orig,
					cutoff_est,
					bw_est,
					K,
					3,
					islog,
					zbar_est
				)

				Fdelta =
					d_bmodel_row23_ddelta(
						delta,
						cutoff_orig,
						bw_orig,
						cutoff_est,
						bw_est,
						K,
						3,
						islog,
						zbar_est
					) * beta'

				if (abs(Fdelta) < 1e-12 | Fdelta >= .) {
					st_numscalar("r(pb_failcode)", 303)
					return
				}

				ddelta_dbeta = -Rbmod / Fdelta

				ddelta_dtheta =
					ddelta_dbeta,
					J(1, Kb, 0),
					1/Fdelta

				Gq = J(Kb, 2*Kb + 1, 0)

				Gq[., 1..Kb] = (
					-hmap.dgamma_dbeta
					- hmap.dgamma_ddelta * ddelta_dtheta[1, 1..Kb]
				)

				Gq[., (Kb+1)..(2*Kb)] = I(Kb)

				Gq[., 2*Kb + 1] =
					-hmap.dgamma_ddelta * ddelta_dtheta[1, 2*Kb + 1]
			}
			else {
				_error(3498, "polbunch_wald_from_unrestricted only handles estimators 1, 2, and 3")
			}

			/*
				Wald statistic.
			*/
			Vq = Gq * V * Gq'

			W = q' * pinv(Vq) * q
			df = rows(q)
			pval = chi2tail(df, W)

			st_numscalar("r(pb_wald)", W)
			st_numscalar("r(pb_p)", pval)
			st_numscalar("r(pb_df)", df)
			st_numscalar("r(pb_delta_U)", delta)
		}

real rowvector saez_transform(
    real rowvector theta,
    real scalar zstarorig,
    real scalar bworig,
    real scalar t0,
    real scalar t1,
    real scalar islog,
    real scalar hastax,
    real scalar dograd,
    real scalar useconstant,
    string scalar Gname
)
{
    real scalar hminus, hplus, B, s, taxratio, L
    real scalar excess_mass, shift, marginal_response, elasticity
    real scalar A, q, disc, x, dlogz, Fx
    real rowvector out, dshift, dmr, dx
    real matrix G

    out = J(1, 6 + hastax, .)
    if (dograd) st_matrix(Gname, J(6 + hastax, 3, .))

    if (cols(theta) < 3) return(out)

    hminus = theta[1]
    hplus  = theta[2]
    B      = theta[3]

    s = hminus + hplus
    if (hminus <= 0 | hplus <= 0 | s <= 0) return(out)
    if (zstarorig <= 0 | bworig <= 0) return(out)

    // excess mass is in bins
    excess_mass = 2 * B / s

    if (islog) {
        if (useconstant) {
            /*
                Constant-density approximation: treat h0:_cons (hminus)
                alone as the counterfactual reference, exactly as the
                other estimators use hminus/h0(z*) alone -- Bsaez itself
                (built from both hminus and hplus) is untouched; only the
                B -> response step drops hplus. At hminus==hplus this
                coincides exactly with the exact log formula below (no
                approximation error at all, unlike the level case),
                since dlogz is already linear in B with no further
                nonlinearity to solve.
            */
            dlogz = B * bworig / hminus
            x     = exp(dlogz)

            shift             = x - 1
            marginal_response = dlogz

            if (dograd) {
                dmr = (-B*bworig/hminus^2, 0, bworig/hminus)
                dshift = x * dmr
            }
        }
        else {
            // hminus/hplus/B are bin counts, so convert bin-width mass to log distance
            dlogz = 2 * B * bworig / s
            x     = exp(dlogz)

            shift             = x - 1
            marginal_response = dlogz

            if (dograd) {
                dmr = (-2*B*bworig/s^2, -2*B*bworig/s^2, 2*bworig/s)
                dshift = x * dmr
            }
        }
    }
    else {
        if (useconstant) {
            // Same idea in levels: r_const = B*bworig/hminus, dropping
            // hplus (and hence the quadratic in x) entirely -- coincides
            // with the exact solve below only to first order in the
            // response, even at hminus==hplus, since the level relation
            // is genuinely nonlinear in x (unlike logs).
            A     = 2 * B * bworig / zstarorig
            shift = A / (2*hminus)
            x     = 1 + shift
            marginal_response = zstarorig * shift

            if (dograd) {
                dshift = (-A/(2*hminus^2), 0, bworig/(hminus*zstarorig))
                dmr    = zstarorig * dshift
            }
        }
        else {
            // B = zstarorig/(2*bworig) * (x-1)*(hminus + hplus/x)
            A    = 2 * B * bworig / zstarorig
            q    = hplus - hminus - A
            disc = q^2 + 4*hminus*hplus
            if (disc < 0) return(out)

            x = (-q + sqrt(disc)) / (2*hminus)
            if (x <= 0) return(out)

            shift             = x - 1
            marginal_response = zstarorig * shift

            if (dograd) {
                Fx = zstarorig/(2*bworig) * (hminus + hplus/x^2)
                if (Fx == 0) return(out)

                dx = J(1,3,0)
                dx[1] = -(zstarorig/(2*bworig)) * (x-1) / Fx
                dx[2] = -(zstarorig/(2*bworig)) * (x-1) / x / Fx
                dx[3] =  1 / Fx

                dshift = dx
                dmr    = zstarorig * dshift
            }
        }
    }

    if (hastax) {
        taxratio = (1-t0)/(1-t1)
        elasticity = .

        if (taxratio > 0 & taxratio != 1) {
            L = ln(taxratio)
            if (islog) elasticity = marginal_response / L
            else       elasticity = ln(x) / L
        }

        out = (hminus, hplus, B, excess_mass, shift, marginal_response, elasticity)
    }
    else {
        out = (hminus, hplus, B, excess_mass, shift, marginal_response)
    }

    if (dograd) {
        G = J(6 + hastax, 3, 0)

        G[1,1] = 1
        G[2,2] = 1
        G[3,3] = 1

        G[4,1] = -2 * B / s^2
        G[4,2] = -2 * B / s^2
        G[4,3] =  2 / s

        G[5,.] = dshift
        G[6,.] = dmr

        if (hastax) {
            if (taxratio <= 0 | taxratio == 1) {
                G[7,.] = J(1,3,.)
            }
            else {
                if (islog) G[7,.] = dmr / L
                else       G[7,.] = dshift / (x * L)
            }
        }

        st_matrix(Gname, G)
    }

    return(out)
}



		end