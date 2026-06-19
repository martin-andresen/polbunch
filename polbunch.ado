				*! polbunch version date 20260618
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
					nodrop ///
					notransform ///
					positive ///
					nonormalize ///
					BOOTreps(integer 500) ///
					vce(string) ///
					log ///
					constant ///
					nodots /// suppress dots for bootstrap progress
					test(string) ///
					nosmallsample ///
					nobayes ///
					nozero ///
					saveunres(string) ///
					Bmodel ///
					wald ///
					norankred ///
					]
					
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
						
						if "`test'"=="" {
							if inlist(`estimator',1,4) loc test wald
							else loc test hausman
						}
						else {
							if !inlist("`test'","none","minimumdistance","wald","hausman","all") {
								noi di as error "Test() can only take wald, minimumdistance, hausman,"all" or none."
								exit 301
							}
							if "`test'"!="wald" & inlist(`estimator',1,4) {
								noi di as error "Only test(wald) supported for estimator 1 and 4 - simple linear restrictions."
								exit 301
							}
						}
						
						

						if `estimator'==4 loc polynomial=0

						tempvar touse
						marksample touse
						preserve
						drop if !`touse'
						
						if "`vce'"=="" loc vce analytic
						else {
							if !inlist("`vce'","analytic","bootstrap","bayes","none") {
								noi di as error "vce() can only contain none, analytic, bootstrap or bayes."
								exit 301
							}
						}
						local grad = cond("`vce'"=="analytic","","nograd")
						
						if `bootreps'<=1 {
							noi di as error "Option bootreps can only take an integer >1 (binned bootstrap)."
							exit 301
						}
						if `polynomial'<0 {
							noi di as error "Polynomial must be a nonnegative integer"
							exit 301
						}
						
						
						// check varlist vs bw opts
						loc nvars: word count `varlist' 
						if (`nvars'==1&"`bw'"=="")|(`nvars'==2&"`bw'"!="") {
							noi di as error "Varlist must either contain 1 variable (earnings z) and option bw be specified (individual level data) or 2 variables (frequency y and earnings bin z) and option bw not be specified (pre-binned data)."
							exit 301
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
										
						//gen dummies
						tempvar fw dupe dum dum2 cons
						gen byte `dum' = `z' > `cutoff_est'
						gen byte `dum2' = `dum'
						gen `cons'=1
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
											local nmiss = 0
										}
									}
								}
							}
								if "`note'"=="note" {
									noi di as text "Note: Polynomial order lowered to `polynomial' because of multicollinearity problems with the specified polynomial."
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
						if inlist("`vce'","bootstrap","bayes") {
							tempname p yorig
							if "`vce'"=="bayes" gen double `p'=`y'/`N'
							else {
								gen double `yorig'=`y'
								recast double `y'
							}
						}
								
						local dotest = inlist(`estimator', 1, 2, 3,4) & "`test'" != "none" & "`vce'"!="none"
						
						//ESTIMATION AND INFERENCE
						tempname b V bs bb VV bmain Vmain b0 V0 b0s bR_raw GR_raw y_raw bU_raw GU_raw d_raw ds Dmain VD
						
						if inlist("`vce'","none","analytic") loc stop=0
						else loc stop=`bootreps'
						forvalues s=0/`stop' {
							if `s'==1&"`dots'"!="nodots" nois _dots 0, title("Performing bootstrap repetitions...") reps(`bootreps')
							if `s'>0 { //resample outcome
								if "`vce'"=="bayes" {
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
								else {
									replace `y'=rgamma(`yorig',1) 
									su `y'
									replace `y'=`y'*`N'/r(sum)
								}
								
							}
							
							//estimate model: single stacked profile branch for estimators 0/1/2/3
						
								if `estimator' == 4 {
									bunch_saez `y' `z' `side' `bunch', ///
										cutoff_orig(`cutoff_orig') ///
										bw_orig(`bw_orig') ///
										zl_excl_orig(`zL_excl_orig') ///
										zh_excl_orig(`zH_excl_orig') ///
										vce(`vce')
								}
								else {
									bunch_profile `y' `z' `side' `bunch', ///
										estimator(`estimator') k(`polynomial') ///
										cutoff_orig(`cutoff_orig') bw_orig(`bw_orig') ///
										cutoff_est(`cutoff_est') bw_est(`bw_est') ///
										l(`L') h(`H') ///
										`log' `normalize' vce(`vce') ///
										initdelta(`dstart') ///
										zbar_est(`zbar_est') ///
										zl_excl_orig(`zL_excl_orig') ///
										zh_excl_orig(`zH_excl_orig') ///
										zl_excl_est(`zL_excl_est') ///
										zh_excl_est(`zH_excl_est') ///
										`positive'
								}

							//STORE RAW RESTRICTED ESTIMATES if using hausman test.
							if `dotest' & "`test'"=="hausman" {
								matrix `bR_raw' = e(b)

								if "`vce'"=="analytic" {
									matrix `GR_raw' = e(G_stack)
									matrix `y_raw'  = e(y_stack)
								}
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
										`bmodel'			
								}
								else {
									saez_transform, zstarorig(`cutoff_orig') bworig(`bw_orig') t0(`t0') t1(`t1') `log' `grad'
									}
								}
								
								matrix `b' = e(b)
								
								if `s' == 0 {
									matrix `bmain' = `b'
									if "`vce'"=="analytic" matrix `Vmain' = e(V)
								}
								else if `bootreps'>1 {
									mat `bs'=nullmat(`bs') \ `b'
								}
								
								
								//IF TESTING: ALSO ESTIMATE UNRESTRICTED MODEL
								if `dotest'&`estimator'!=4 {
									bunch_profile `y' `z' `side' `bunch', ///
										estimator(0) k(`polynomial') ///
										cutoff_orig(`cutoff_orig') bw_orig(`bw_orig') ///
										cutoff_est(`cutoff_est') bw_est(`bw_est') ///
										l(`L') h(`H') ///
										`log' `normalize' vce(`vce') ///
										initdelta(`dstart') ///
										zbar_est(`zbar_est') ///
										zl_excl_orig(`zL_excl_orig') ///
										zh_excl_orig(`zH_excl_orig') ///
										zl_excl_est(`zL_excl_est') ///
										zh_excl_est(`zH_excl_est') ///
										`positive'
									
									if "`test'"=="hausman" {
										matrix `bU_raw' = e(b)

										if "`vce'"=="analytic" {
											matrix `GU_raw' = e(G_stack)
										}

										polbunch_modeldiff, ///
											estimator(`estimator') ///
											k(`polynomial') ///
											bu(`bU_raw') ///
											br(`bR_raw') ///
											cutofforig(`cutoff_orig') ///
											bworig(`bw_orig') ///
											cutoffest(`cutoff_est') ///
											bwest(`bw_est') ///
											zbar(`zbar_est') ///
											`normalize' ///
											`log'

										matrix `d_raw' = r(d)

										if `s' == 0 {
											matrix `Dmain' = `d_raw'
										}
										else if `bootreps' > 1 {
											matrix `ds' = nullmat(`ds') \ `d_raw'
										}
									}

									if `s' == 0 {
										matrix `b0' = e(b)
										if "`vce'"=="analytic" matrix `V0' = e(V)
									}
									else if `bootreps' > 1 {
										matrix `b0s' = nullmat(`b0s') \ e(b)
									}
								}
								if `s' > 0 noi _dots `s' 0
								
							}
							
								
							//bootstrap inference sunmmary & test
							if inlist("`vce'","bootstrap","bayes") {
								clear
								svmat double `bs'
								corr _all, cov
								mat `Vmain'=r(C)
								if `dotest' {
									if "`test'"=="hausman" {
										clear
										svmat double `ds'
										corr _all, cov
										matrix `VD' = r(C)
									}
									else if `estimator'!=4 {
										clear
										svmat `b0s'
										corr _all, cov
										matrix `V0' = r(C)
									}

								}
							}
							

								
							//TEST RESTRICTIONS
			
							if `dotest' {
								if `estimator'==4 { //saez: Post main model to b0 V0
										mat `b0' = `bmain'
										mat `V0' = `Vmain'

									}
								local testname=cond("`test'"=="wald","Wald test",cond("`test'"=="minimumdistance","Minumum-distance test","Hausman test"))
								
								if inlist("`test'","wald","minimumdistance") {
										local nm: colnames `b0'
										local neq: coleq `b0'
										mat colnames `V0'=`nm'
										mat rownames `V0'=`nm'
										mat coleq `V0'=`neq'
										mat roweq `V0'=`neq'
										
										ereturn post `b0' `V0'
										ereturn local properties "b V"
									}
									
								
								if `estimator'!=4 {
									if "`test'"=="minimumdistance" {
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
										
										local delta_md=r(delta)
									}
									else if "`test'"=="wald" {
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
									else if "`test'"=="hausman" {
										if "`vce'"=="analytic" {
											capture noisily polbunch_modeltest, ///
												estimator(`estimator') ///
												k(`polynomial') ///
												bu(`bU_raw') gu(`GU_raw') ///
												br(`bR_raw') gr(`GR_raw') ///
												ystack(`y_raw') ///
												cutofforig(`cutoff_orig') ///
												bworig(`bw_orig') ///
												cutoffest(`cutoff_est') ///
												bwest(`bw_est') ///
												zbar(`zbar_est') ///
												`normalize' ///
												`log'
										}
										else {
											capture noisily polbunch_diff_test, d(`Dmain') v(`VD')
										}
									}
								}
								else { //Saez: Simple Wald test of h0 vs h1
									if "`log'"=="" capture test _b[h1:_cons] - _b[h0:_cons] - (`bw_orig'/`cutoff_orig')*_b[bunching:number_bunchers] = 0
									else capture test _b[h1:_cons] = _b[h0:_cons]
								}

								local test_rc = _rc

								local failcode = 0
								capture confirm scalar r(failcode)
								if !_rc {
									local failcode = r(failcode)
									if missing(`failcode') local failcode = 0
								}

								if (`test_rc' != 0) | (`failcode' != 0) {
									local dotest = 0
									local test_failed = 1

									if `test_rc' != 0 local test_failcode = `test_rc'
									else              local test_failcode = `failcode'

									noi di as text  "Note: model-assumption test could not be computed; `testname' statistic is not reported."
									noi di as text 	"Consider one of the alternatives in test(). `testname' return code = " as result `test_failcode'
								}
								else {
									local chi2  = r(chi2)
									local p_mod = r(p)
									local df    = r(df)
								}
							}

						
						//POST RESULTS
						su `z' if `side'==-1, meanonly
						loc dL=r(mean)
						noi di "`dL'"
						su `z' if `side'==1, meanonly
						loc dR=r(mean)
						noi di "`dR'"
						
						restore
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
							estadd scalar chi2=`chi2'
							estadd scalar p_mod=`p_mod'
							estadd scalar df_mod=`df'
							if "`test'"=="minimumdistance" estadd scalar delta_md = `delta_md'
						}
						ereturn scalar polynomial=`polynomial'
						ereturn scalar bandwidth=`bw'
						ereturn scalar cutoff=`cutoff'
						ereturn scalar lower_limit=`zL_excl_orig'
						ereturn scalar upper_limit=`zH_excl_orig'
						ereturn local normalize="`normalize'"
						ereturn scalar estimator=`estimator'
						if "`vce'"!="none" ereturn local cmd "polbunch"
						ereturn local cmdname "polbunch"
						ereturn local title 	"Polynomial bunching estimates"
						ereturn local cmdline 	"polbunch `0'"
						ereturn matrix table=`table'
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

						if "`vce'"=="analytic"&"`transform'"=="notransform" estadd local vcetype "analytic"
						if "`vce'"=="analytic"&"`transform'"=="" estadd local vcetype "delta method"
						else if "`vce'"=="bootstrap" estadd local vcetype "binned bootstrap"
						else if "`vce'"=="bayes" estadd local vcetype "bayesian bootstrap"
						if "`log'"=="log" ereturn scalar log=1
						else ereturn scalar log=0
			
						//Display results
						noi {
							di _newline
							di "`e(title)'"
							eret di
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
								di as txt "`testname' of model assumptions:" ///
									_col(`=`W'-35') "Chi2(`df') test statistic" ///
									_col(`=`W'-10') as res %10.4f `chi2'
								di as txt _col(`=`W'-35') as txt "p-value" ///
									_col(`=`W'-10') as res %10.4f `p_mod'
								di as txt "{hline `W'}"
							}
							if inlist(`estimator',1,2) {
								di "Note: Estimator is not consistent with iso-elastic labor supply model and thus biased."
							}
							if "`constant'"!=""&"`t0'"!=""&"`t1'"!=""&"`transform'"!="notransform" {
								di "Note: Using the constant approximation to the density to calculate the elasticity may lead to bias."
							}
						}
					
					if "`vce'"!="none" set coeftabresults `coeftabresults'
					}
						
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
				[ nonormalize LOG vce(string) initdelta(real 0.05) positive]

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
				`normalized0', ///
				`islog0', ///
				`zl_excl_orig', ///
				`zh_excl_orig', ///
				`zl_excl_est', ///
				`zh_excl_est', ///
				`zbar_est', ///
				`dovar0', ///
				`initdelta', ///
				`positive0' ///
			)

			matrix `b' = r_b_profile
			matrix `V' = r_V_profile

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
				[ LOG CONSTANT T0(numlist min=1 max=1) T1(numlist min=1 max=1) nograd nonormalize ZBAR(real 0) MASSOBS(real 0) bmodel ]

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
			if "`grad'" != "nograd" {
				capture confirm matrix e(V)
				if !_rc {
					matrix `Vtheta' = e(V)
					local dograd = 1
				}
			}

			/* Flags */
			local islog       = ("`log'"      != "")
			local constant0   = ("`constant'" != "")
			local hastax0 = ("`t0'" != "" & "`t1'" != "")
			local normalized0 = ("`normalize'" != "nonormalize")

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

			/* Call Mata transform */
			local Btype2 = 1
			if "`bmodel'" != "" local Btype2 = 0

			mata: bunch_transform( ///
				st_matrix("`theta'"), ///
				`estimator', ///
				`k', ///
				`cutofforig', ///
				`cutoffest', ///
				`bworig', ///
				`bwest', ///
				`xscale', ///
				`low', ///
				`high', ///
				`islog', ///
				`constant0', ///
				`hastax0', ///
				`t0', ///
				`t1', ///
				`dograd', ///
				`normalized0', ///
				`zbar', ///
				`massobs', ///
				`Btype2', ///
				`zlexcl', ///
				`zhexcl', ///
				`zlexclest', ///
				`zhexclest' ///
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

			if (`estimator' == 2 | (`estimator' == 3 & `constant0')) {
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

			if (`estimator' == 2 | (`estimator' == 3 & `constant0')) {
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
        [T0(numlist max=1) T1(numlist max=1) log nograd]

    tempname b0 V0 theta Vtheta b G V

    matrix `b0' = e(b)
    if (colsof(`b0') < 3) {
        di as err "e(b) must contain theta=(h0:_cons,h1:_cons,B) in columns 1..3"
        exit 503
    }
	
	noi di "`dograd'"
    local islog  = ("`log'" != "")
    local dograd = ("`grad'" != "nograd")

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
        `islog', `hastax', `dograd', "`G'" ///
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
			[vce(string)]

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
		

		mata: saez_run("`y_t'", "`side_t'", "`bunch_t'", `a0', `a1', `dovar0')

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
		}

		ereturn local cmd "bunch_saez"
		ereturn scalar estimator = 4
		ereturn scalar cutoff_orig = `cutoff_orig'
		ereturn scalar bw_orig = `bw_orig'
		ereturn scalar zL_excl_orig = `zl_excl_orig'
		ereturn scalar zH_excl_orig = `zh_excl_orig'
		ereturn scalar saez_a0 = `a0'
		ereturn scalar saez_a1 = `a1'
		ereturn scalar saez_width_excl = `width_excl'
	end

	cap program drop polbunch_modeldiff
	program define polbunch_modeldiff, rclass
		syntax , ///
			ESTimator(integer) ///
			K(integer) ///
			BU(name) ///
			BR(name) ///
			cutofforig(real) ///
			bworig(real) ///
			cutoffest(real) ///
			bwest(real) ///
			zbar(real) ///
			[ nonormalize log ]

		local normalized0 = ("`normalize'" != "nonormalize")
		local islog0      = ("`log'" != "")

		mata: polbunch_modeldiff_mata( ///
			"`bu'", "`br'", ///
			`estimator', `k', ///
			`cutofforig', `bworig', `cutoffest',`bwest',`zbar', ///
			`normalized0', `islog0' ///
		)

		matrix d = r(pb_model_d)
		return matrix d = d
		return scalar failcode = r(pb_modeldiff_failcode)
	end

	cap program drop polbunch_diff_test
	program define polbunch_diff_test, rclass
		syntax , D(name) V(name)

		mata: polbunch_diff_test_mata("`d'", "`v'")

		return scalar chi2 = r(pb_diff_chi2)
		return scalar p    = r(pb_diff_p)
		return scalar df   = r(pb_diff_df)
		return scalar failcode = r(pb_diff_failcode)
	end

	cap program drop polbunch_modeltest
	program define polbunch_modeltest, rclass
		version 16.0

		syntax , ///
			ESTimator(integer) ///
			K(integer) ///
			BU(name) GU(name) ///
			BR(name) GR(name) ///
			YSTACK(name) ///
			cutofforig(real) ///
			bworig(real) ///
			cutoffest(real) ///
			bwest(real) ///
			ZBAR(real) ///
			[ NONORMALIZE LOG ]

		local normalized0 = ("`normalize'" != "nonormalize")
		local islog0      = ("`log'" != "")

		mata: polbunch_modeltest_mata( ///
			"`bu'", "`gu'", ///
			"`br'", "`gr'", ///
			"`ystack'", ///
			`estimator', ///
			`k', ///
			`cutofforig', ///
			`bworig', ///
			`cutoffest', ///
			`bwest', ///
			`zbar', ///
			`normalized0', ///
			`islog0' ///
		)

		return scalar chi2 = r(pb_model_chi2)
		return scalar p    = r(pb_model_p)
		return scalar df   = r(pb_model_df)
		return scalar failcode = r(pb_model_failcode)
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

			local normalized0 = ("`normalize'" != "nonormalize")
			local islog0      = ("`log'" != "")

			mata: polbunch_wald_from_unrestricted( ///
				"`b'", ///
				"`V'", ///
				`estimator', ///
				`cutofforig', ///
				`bworig', ///
				`cutoffest', ///
				`bwest', ///
				`k', ///
				`normalized0', ///
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

		local normalized0 = ("`normalize'" != "nonormalize")
		local islog0      = ("`log'" != "")
		local positive0   = ("`positive'" != "")

		capture noisily mata: polbunch_mdt_mata( ///
			"`b'", ///
			"`V'", ///
			`estimator', ///
			`cutofforig', ///
			`bworig', ///
			`cutoffest', ///
			`bwest', ///
			`k', ///
			`normalized0', ///
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
		real scalar cutoff_est,
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
		real scalar cutoff_est,
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
		real scalar dovar
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
			Vout = variance_multinomial(X, ystack, 0)
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
		real scalar normalized,
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
		real scalar normalized,
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
			normalized,
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
				normalized,
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
		real scalar normalized,
		real scalar islog,
		real scalar dograd,
		real scalar zL_excl_est,
		real scalar zH_excl_est
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
			normalized,
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
				normalized,
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
			real scalar normalized,
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
				/* Chetty restriction: bw * B = delta * int_{zstar}^{zbar} h0(z) dz */
				R = (delta / bw_est) * intbasis(cutoff_est, zbar_est, K)
			}
			else if (estimator == 3) {
				/* Theoretically consistent restriction: bw * B = int_{zstar}^{zstar+r(delta)} h0(z) dz */
				r = response_length(delta, cutoff_orig, bw_orig, cutoff_est,bw_est,islog)
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
		real scalar normalized,
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
				normalized,
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
		real scalar normalized,
		real scalar islog,
		real scalar zbar_est
	)
	{
		real scalar r, dr

		if (1 + delta <= 0) {
			_error(3498, "delta must be greater than -1")
		}

		if (estimator == 2) {
			return(intbasis(cutoff_est, zbar_est, K) / bw_est)
		}

		if (estimator == 3) {
			r  = response_length(delta, cutoff_orig, bw_orig,cutoff_est,bw_est, islog)
			dr = d_response_length_ddelta(delta, cutoff_orig, bw_orig,cutoff_est,bw_est,islog)

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
			real scalar normalized,
			real scalar islog,
			real scalar zL_excl_orig,
			real scalar zH_excl_orig,
			real scalar zL_excl_est,
			real scalar zH_excl_est,
			real scalar ntheta
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
				Rlo = intbasis(zL_excl_est, cutoff_est, K) / bw_est
				Rhi = intbasis(cutoff_est, zH_excl_est, K) / ((1 + delta) * bw_est)
				R[1, 1..Kb] = Rlo + Rhi
			}
			else if (estimator == 3) {
				Rlo = intbasis(zL_excl_est, cutoff_est, K) / bw_est
				h1map = h1coef_map(
					J(1, K+1, 0),
					delta,
					estimator,
					K,
					cutoff_orig,
					bw_orig,
					cutoff_est,
					bw_est,
					normalized,
					islog,
					1
				)
				Rhi = (intbasis(cutoff_est, zH_excl_est, K) * h1map.dgamma_dbeta) / bw_est
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
		real scalar normalized,
		real scalar islog,
		real scalar zL_excl_orig,
		real scalar zH_excl_orig,
		real scalar zL_excl_est,
		real scalar zH_excl_est,
		real scalar ntheta
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
				normalized,
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
			real scalar normalized,
			real scalar islog,
			real scalar zbar_est,
			real scalar ntheta,
			real scalar zL_excl_est,
			real scalar zH_excl_est
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
				R[1, 1..(K+1)] = (delta / bw_est) * intbasis(cutoff_est, zbar_est, K)
			}
			else if (estimator == 3) {
				r = response_length(delta, cutoff_orig, bw_orig, cutoff_est,bw_est,islog)
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
			real scalar normalized,
			real scalar islog,
			real scalar zL_excl_orig,
			real scalar zH_excl_orig,
			real scalar zL_excl_est,
			real scalar zH_excl_est,
			real scalar zbar_est,
			real scalar dograd
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

				Xcf = cf_mass_row(delta, cutoff_orig, bw_orig,cutoff_est,bw_est, K, estimator, normalized, islog, zL_excl_orig, zH_excl_orig,zL_excl_est,zH_excl_est, ntheta)
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

				Xcf = cf_mass_row(delta, cutoff_orig, bw_orig,cutoff_est,bw_est, K, estimator, normalized, islog, zL_excl_orig, zH_excl_orig,zL_excl_est,zH_excl_est, ntheta)
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
				h1 = h1design23(delta, zR, cutoff_orig, bw_orig,cutoff_est,bw_est, K, estimator, normalized, islog, dograd,zL_excl_est,zH_excl_est)

				ntheta = Kb
				Xcf    = cf_mass_row(delta, cutoff_orig, bw_orig,cutoff_est,bw_est, K, estimator, normalized, islog, zL_excl_orig, zH_excl_orig,zL_excl_est,zH_excl_est, ntheta)
				Xbmod  = bmodel_row(delta, cutoff_orig, bw_orig,cutoff_est,bw_est, K, estimator, normalized, islog, zbar_est, ntheta,zL_excl_est,zH_excl_est)
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
						normalized,
						islog,
						zL_excl_orig,
						zH_excl_orig,
						zL_excl_est,
						zH_excl_est,
						ntheta
					)

					dXbmod = d_bmodel_row_ddelta(
						delta,
						cutoff_orig,
						bw_orig,
						cutoff_est,
						bw_est,
						K,
						estimator,
						normalized,
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
			real scalar normalized,
			real scalar islog,
			real scalar zL_excl_orig,
			real scalar zH_excl_orig,
			real scalar zL_excl_est,
			real scalar zH_excl_est,
			real scalar zbar_est
		)
		{
			real colvector ystack, theta
			struct design_out scalar D

			ystack = make_ystack(y, side, bunch, estimator, Hstar_obs)

			D = make_design(delta, z, side, bunch, cutoff_orig, bw_orig, cutoff_est,bw_est,K, estimator, normalized, islog, zL_excl_orig, zH_excl_orig,zL_excl_est,zH_excl_est, zbar_est, 0)

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
			real scalar normalized,
			real scalar islog,
			real scalar zL_excl_orig,
			real scalar zH_excl_orig,
			real scalar zL_excl_est,
			real scalar zH_excl_est,
			real scalar zbar_est,
			real scalar positive
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
				normalized,
				islog,
				zL_excl_orig,
				zH_excl_orig,
				zL_excl_est,
				zH_excl_est,
				zbar_est,
				0
			)

			theta = qrsolve(D.X, ystack)
			resid = ystack - D.X * theta

			return(quadcross(resid, resid))
		}


		// -----------------------------------------------------------------------------
		// Variance estimation
		// -----------------------------------------------------------------------------


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
			real scalar L,
			real scalar H,
			real scalar islog,
			real scalar constant,
			real scalar hastax,
			real scalar t0,
			real scalar t1,
			real scalar dograd,
			real scalar normalized,
			real scalar zbar_est,
			real scalar Hstar_obs,
			real scalar Btype2,
			real scalar zL_excl_orig,
			real scalar zH_excl_orig,
			real scalar zL_excl_est,
			real scalar zH_excl_est
		)
		{
			struct hcoef_out scalar h1

			real rowvector xcut, dm, dB, RB, dRB
			real rowvector beta, gamma
			real rowvector dr, dF_dbeta, dMR, dshift
			real matrix G
			real rowvector b

			real scalar Kb, nout, i, hasresp
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

			nout =
				2*Kb +                                      // h0,h1
				2 +                                         // B, EM
				(estimator==2 | (estimator==3 & constant)) + // delta
				2 +                                         // shift, MR
				hastax

			b = J(1, nout, .)
			if (dograd) G = J(nout, cols(theta), 0)

			// output indices
			obeta  = ibeta
			ogamma = (Kb+1)..(2*Kb)
			oB     = 2*Kb + 1
			oEM    = 2*Kb + 2

			if (estimator==2 | (estimator==3 & constant)) {
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
					normalized,
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
			if (estimator == 0 | estimator == 1) {
				b[1,oB] = B
				if (dograd) G[oB,iBraw] = 1
			}
			else if (estimator == 2) {
				if (Btype2 == 0) {
					/* Estimator 2, model-implied B: B_model = delta/bw * int_{z*}^{zbar} h0(z) dz */
					RB = bmodel_row23(delta, cutoff_orig, bw_orig,cutoff_est,bw_est, K, estimator, normalized, islog, zbar_est)

					B = RB * beta'
					b[1,oB] = B

					if (dograd) {
						G[oB, ibeta] = RB
						dRB = d_bmodel_row23_ddelta(delta, cutoff_orig, bw_orig,cutoff_est,bw_est, K, estimator, normalized, islog, zbar_est)
						G[oB, idelta] = dRB * beta'
					}
				}
				else if (Btype2 == 1) {
					/* Estimator 2, Chetty-style reduced-form B: B_reduced = Hstar_obs - int_{lower}^{upper} h0(z) dz / bw */
					RB = intbasis(zL_excl_est,zH_excl_est, K) / bw_est
					B = Hstar_obs - RB * beta'
					b[1,oB] = B

					if (dograd) {
						G[oB, ibeta]  = -RB
						G[oB, idelta] = 0
					}
				}
				else {
					_error(3498, "Btype2 must be 0 for B_model or 1 for B_reduced")
				}
			}
			else if (estimator == 3) {
				/* Estimator 3 always reports the theoretically consistent model-implied B. */
				RB = bmodel_row23(delta, cutoff_orig, bw_orig,cutoff_est, bw_est, K, estimator, normalized, islog, zbar_est)

				B = RB * beta'
				b[1,oB] = B

				if (dograd) {
					G[oB, ibeta] = RB
					dRB = d_bmodel_row23_ddelta(delta, cutoff_orig, bw_orig,cutoff_est,bw_est, K, estimator, normalized, islog, zbar_est)
					G[oB, idelta] = dRB * beta'
				}
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

			// delta for estimator 2 and estimator 3 with constant approximation
			if (estimator == 2 | (estimator == 3 & constant)) {
				b[1,odelta] = delta
				if (dograd) G[odelta,idelta] = 1
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
			else if (estimator == 3) {
				// estimator 3: shift is structural delta
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
				// estimators 0/1/2: solve integral equation using eresp()
				r = eresp(B, cutoff_est, beta, bw_est, xscale)

				if (r >= .) {
					/*
						If no admissible root exists for the missing-mass equation,
						report coefficients, B, EM, and delta if present, but drop
						shift/marginal_response/elasticity from the posted vector.
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
		void profQ_opt(
			real scalar todo,
			real rowvector p,
			real colvector y,
			real colvector z,
			real colvector side,
			real colvector bunch,
			real rowvector pars,
			real scalar val,
			real rowvector grad,
			real matrix hess
		)
		{
			real scalar Hstar_obs, cutoff_orig, bw_orig,cutoff_est,bw_est, K
			real scalar estimator, normalized, islog
			real scalar zL_excl_orig, zH_excl_orig, zbar_est, positive
			real scalar delta

			Hstar_obs    = pars[1]
			cutoff_orig  = pars[2]
			bw_orig      = pars[3]
			cutoff_est 	 = pars[4]
			bw_est 		 = pars[5]
			K            = pars[6]
			estimator    = pars[7]
			normalized   = pars[8]
			islog        = pars[9]
			zL_excl_orig = pars[10]
			zH_excl_orig = pars[11]
			zL_excl_est = pars[12]
			zH_excl_est = pars[13]
			zbar_est     = pars[14]
			positive     = pars[15]

			delta = p[1]

			val = profQ(
				delta,
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
				normalized,
				islog,
				zL_excl_orig,
				zH_excl_orig,
				zL_excl_est,
				zH_excl_est,
				zbar_est,
				positive
			)
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
			real scalar normalized,
			real scalar islog,
			real scalar zL_excl_orig,
			real scalar zH_excl_orig,
			real scalar zL_excl_est,
			real scalar zH_excl_est,
			real scalar zbar_est,
			real scalar dovar,
			real scalar initdelta,
			real scalar positive

		)
		{ 
			real colvector y, z, side, bunch
			real scalar Kb, Hstar_obs, lndelta_hat, delta_hat, gi, lb
			real rowvector theta_hat, beta_hat, b, pars, phat, dgrid, qgrid
			real matrix Vout
			transmorphic S
			
			struct design_out scalar D
			real colvector ystack, mu
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
					normalized,
					islog,
					zL_excl_orig,
					zH_excl_orig,
					zL_excl_est,
					zH_excl_est,
					zbar_est
				)
				b = theta_hat
			}
			else if (estimator == 2 | estimator == 3) {
				pars = (
					Hstar_obs,
					cutoff_orig,
					bw_orig,
					cutoff_est,
					bw_est,
					K,
					estimator,
					normalized,
					islog,
					zL_excl_orig,
					zH_excl_orig,
					zL_excl_est,
					zH_excl_est,
					zbar_est,
					positive
				)

				S = optimize_init()
				optimize_init_evaluator(S, &profQ_opt())
				optimize_init_evaluatortype(S, "d0")
				optimize_init_which(S, "min")
				optimize_init_conv_maxiter(S, 200)

				if (positive == 1) {
					lb = 1e-8
					if (initdelta <= lb | initdelta >= .) initdelta = 0.05
				}
				else {
					lb = -1 + 1e-8
					if (initdelta <= lb | initdelta >= .) initdelta = 0.05
				}

				optimize_init_params(S, initdelta)

				optimize_init_argument(S, 1, y)
				optimize_init_argument(S, 2, z)
				optimize_init_argument(S, 3, side)
				optimize_init_argument(S, 4, bunch)
				optimize_init_argument(S, 5, pars)

				phat = optimize(S)
				delta_hat = phat[1]
				
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
					normalized,
					islog,
					zL_excl_orig,
					zH_excl_orig,
					zL_excl_est,
					zH_excl_est,
					zbar_est
				)

				theta_hat = beta_hat, delta_hat
				b = theta_hat
			}
			else {
				_error(3498, "profile_run only handles estimators 0, 1, 2, and 3")
			}

			if (dovar == 1) {
				D = make_design(delta_hat, z, side, bunch,
					cutoff_orig, bw_orig,cutoff_est,bw_est, K, estimator, normalized, islog,
					zL_excl_orig, zH_excl_orig,zL_excl_est,zH_excl_est, zbar_est, 1)

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

				Vout = variance_multinomial(Gv, ystack, 0)
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
		real scalar normalized,
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
			normalized,
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
				normalized,
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
						normalized,
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
		real scalar normalized,
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
			normalized,
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
			normalized,
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
		real scalar normalized,
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
			normalized,
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
		real scalar normalized,
		real scalar islog,
		real scalar zbar_est,
		real scalar positive,
		real scalar initdelta
	)
	{
		real scalar Kb, ntheta, df, pval
		real scalar delta_hat, W_hat, W0
		real scalar lo, hi, center, span
		real scalar i, j, d, W, bestd, bestW
		real scalar a, b, c, x1, x2, f1, f2, gr, iter
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
				normalized,
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

			if (positive == 1) {
				if (d <= lo) continue
			}
			else {
				if (d <= lo) continue
			}

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
				normalized,
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
			normalized,
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
			normalized,
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
					normalized,
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
					normalized,
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

	void polbunch_diff_test_mata(string scalar dname, string scalar Vname)
	{
		real rowvector d
		real matrix V
		real scalar stat, df, pval

		st_numscalar("r(pb_diff_chi2)", .)
		st_numscalar("r(pb_diff_p)", .)
		st_numscalar("r(pb_diff_df)", .)
		st_numscalar("r(pb_diff_failcode)", 0)

		d = st_matrix(dname)
		V = st_matrix(Vname)

		if (missing(d) | missing(V)) {
			st_numscalar("r(pb_diff_failcode)", 101)
			return
		}

		if (rows(d) != 1) d = d'

		if (rows(V) != cols(d) | cols(V) != cols(d)) {
			st_numscalar("r(pb_diff_failcode)", 102)
			return
		}

		df = rank(V)
		if (df <= 0) {
			st_numscalar("r(pb_diff_failcode)", 103)
			return
		}

		stat = d * pinv(V) * d'
		pval = chi2tail(df, stat)

		st_numscalar("r(pb_diff_chi2)", stat)
		st_numscalar("r(pb_diff_p)", pval)
		st_numscalar("r(pb_diff_df)", df)
	}

	void polbunch_modeldiff_mata(
		string scalar bUname,
		string scalar bRname,
		real scalar estimator,
		real scalar K,
		real scalar cutoff_orig,
		real scalar bw_orig,
		real scalar cutoff_est,
		real scalar bw_est,
		real scalar zbar_est,
		real scalar normalized,
		real scalar islog
	)
	{
		real scalar Kb
		real rowvector bU, bR, beta, g, RB
		struct hcoef_out scalar hmap

		st_numscalar("r(pb_modeldiff_failcode)", 0)
		st_matrix("r(pb_model_d)", J(1, 1, .))

		Kb = K + 1
		bU = st_matrix(bUname)
		bR = st_matrix(bRname)

		if (estimator == 2 | estimator == 3) {
			beta = bR[1, 1..Kb]

			hmap = h1coef_map(
				beta,
				bR[1, Kb+1],
				estimator,
				K,
				cutoff_orig,
				bw_orig,
				cutoff_est,
				bw_est,
				normalized,
				islog,
				0
			)

			RB = bmodel_row23(
				bR[1, Kb+1],
				cutoff_orig,
				bw_orig,
				cutoff_est,
				bw_est,
				K,
				estimator,
				normalized,
				islog,
				zbar_est
			)

			g = beta, hmap.gamma, RB * beta'
		}
		else {
			st_numscalar("r(pb_modeldiff_failcode)", 201)
			return
		}

		if (cols(bU) != cols(g)) {
			st_numscalar("r(pb_modeldiff_failcode)", 202)
			return
		}

		st_matrix("r(pb_model_d)", bU - g)
	}

	void polbunch_modeltest_mata(
		string scalar bUname,
		string scalar GUname,
		string scalar bRname,
		string scalar GRname,
		string scalar yname,
		real scalar estimator,
		real scalar K,
		real scalar cutoff_orig,
		real scalar bw_orig,
		real scalar cutoff_est,
		real scalar bw_est,
		real scalar zbar_est,
		real scalar normalized,
		real scalar islog
	)
	{
		real scalar Kb, N, df, stat, pval
		real rowvector bU, bR, beta, g
		real colvector y
		real matrix GU, GR, AU, AR, Jg, Ad, Vm, Vd
		struct hcoef_out scalar hmap
		real rowvector RB, dRB

		st_numscalar("r(pb_model_chi2)", .)
		st_numscalar("r(pb_model_p)", .)
		st_numscalar("r(pb_model_df)", .)
		st_numscalar("r(pb_model_failcode)", 0)

		Kb = K + 1

		bU = st_matrix(bUname)
		GU = st_matrix(GUname)
		bR = st_matrix(bRname)
		GR = st_matrix(GRname)
		y  = st_matrix(yname)

		if (rows(y) == 1) y = y'

		if (missing(bU) | missing(GU) | missing(bR) | missing(GR) | missing(y)) {
			st_numscalar("r(pb_model_failcode)", 101)
			return
		}

		if (rows(GU) != rows(y) | rows(GR) != rows(y)) {
			st_numscalar("r(pb_model_failcode)", 102)
			return
		}

		N = sum(y)
		if (N <= 0 | N >= .) {
			st_numscalar("r(pb_model_failcode)", 103)
			return
		}

		Vm = diag(y) - (y * y') / N
		Vm = (Vm + Vm') / 2

		real rowvector sU, sR
		real matrix GUs, GRs

		sU = sqrt(colsum(GU:^2))
		sR = sqrt(colsum(GR:^2))

		if (any(sU :<= 0) | any(sU :>= .) |
			any(sR :<= 0) | any(sR :>= .)) {
			st_numscalar("r(pb_model_failcode)", 104)
			return
		}

		GUs = GU :/ sU
		GRs = GR :/ sR

		AU = diag(1 :/ sU) * pinv(GUs)
		AR = diag(1 :/ sR) * pinv(GRs)

		if (estimator == 1) {
			/*
				U: (beta, gamma, B)
				R: (beta, B)
				g(beta,B) = (beta, beta, B)
			*/
			if (cols(bU) != 2*Kb + 1 | cols(bR) != Kb + 1) {
				st_numscalar("r(pb_model_failcode)", 201)
				return
			}

			beta = bR[1, 1..Kb]
			g = beta, beta, bR[1, Kb+1]

			Jg = J(2*Kb + 1, Kb + 1, 0)
			Jg[1..Kb, 1..Kb]              = I(Kb)
			Jg[(Kb+1)..(2*Kb), 1..Kb]     = I(Kb)
			Jg[2*Kb+1, Kb+1]              = 1
		}
		else if (estimator == 2 | estimator == 3) {
			/*
				U: (beta, gamma, B)
				R: (beta, delta)
				g(beta,delta) = (beta, gamma(beta,delta), B(beta,delta))
			*/
			if (cols(bU) != 2*Kb + 1 | cols(bR) != Kb + 1) {
				st_numscalar("r(pb_model_failcode)", 202)
				return
			}

			beta = bR[1, 1..Kb]

			hmap = h1coef_map(
				beta,
				bR[1, Kb+1],
				estimator,
				K,
				cutoff_orig,
				bw_orig,
				cutoff_est,
				bw_est,
				normalized,
				islog,
				1
			)

			RB = bmodel_row23(
				bR[1, Kb+1],
				cutoff_orig,
				bw_orig,
				cutoff_est,
				bw_est,
				K,
				estimator,
				normalized,
				islog,
				zbar_est
			)
			dRB = d_bmodel_row23_ddelta(
				bR[1, Kb+1],
				cutoff_orig,
				bw_orig,
				cutoff_est,
				bw_est,
				K,
				estimator,
				normalized,
				islog,
				zbar_est
			)

			g = beta, hmap.gamma, RB * beta'

			Jg = J(2*Kb + 1, Kb + 1, 0)

			Jg[1..Kb, 1..Kb] = I(Kb)

			Jg[(Kb+1)..(2*Kb), 1..Kb] = hmap.dgamma_dbeta
			Jg[(Kb+1)..(2*Kb), Kb+1]  = hmap.dgamma_ddelta

			Jg[2*Kb+1, 1..Kb] = RB
			Jg[2*Kb+1, Kb+1]  = dRB * beta'
		}
		else if (estimator == 4) {
			/*
				Saez:
				U: (h0, h1, B)
				R: (h, B)
				g(h,B) = (h, h, B)

				This requires a restricted Saez estimate with two parameters.
			*/
			if (cols(bU) != 3 | cols(bR) != 2) {
				st_numscalar("r(pb_model_failcode)", 204)
				return
			}

			g = bR[1,1], bR[1,1], bR[1,2]

			Jg = (1, 0 \ 
				  1, 0 \ 
				  0, 1)
		}
		else {
			st_numscalar("r(pb_model_failcode)", 205)
			return
		}

		real matrix C, Vq, Vqs
		real rowvector q, sq, qs
		real scalar qdf

		Ad = AU - Jg * AR
		Vd = Ad * Vm * Ad'
		Vd = (Vd + Vd') / 2

		real matrix Jgs, Uj, Vtj
		real colvector sj
		real rowvector scaleJ
		real scalar rankJ

		/*
			Column scaling improves the SVD numerically but does not alter
			the column space of Jg.
		*/
		scaleJ = sqrt(colsum(Jg:^2))

		if (any(scaleJ :<= 0) | any(scaleJ :>= .)) {
			st_numscalar("r(pb_model_failcode)", 105)
			return
		}

		Jgs = Jg :/ scaleJ

		/*
			Full SVD is needed because the thin SVD would return only
			cols(Jg) left singular vectors and omit the orthogonal complement.
		*/
		fullsvd(Jgs, Uj, sj, Vtj)

		rankJ = rank(Jgs)
		qdf  = rows(Jg) - rankJ

		if (rankJ != cols(Jg) | qdf <= 0) {
			st_numscalar("r(pb_model_failcode)", 105)
			return
		}

		/*
			The last qdf columns of Uj span the left null space:
				C' * Jg = 0.
		*/
		C = Uj[., (rankJ + 1)..rows(Jg)]

		if (cols(C) != qdf) {
			st_numscalar("r(pb_model_failcode)", 105)
			return
		}
		q  = (bU - g) * C
		Vq = C' * Vd * C
		Vq = (Vq + Vq') / 2

		sq = sqrt(diagonal(Vq))'

		if (any(sq :<= 0) | any(sq :>= .)) {
			st_numscalar("r(pb_model_failcode)", 106)
			return
		}

		qs  = q :/ sq
		Vqs = Vq :/ (sq' * sq)
		Vqs = (Vqs + Vqs') / 2

		stat = qs * pinv(Vqs) * qs'
		df   = qdf
		pval = chi2tail(df, stat)

		st_numscalar("r(pb_model_chi2)", stat)
		st_numscalar("r(pb_model_p)", pval)
		st_numscalar("r(pb_model_df)", df)
		st_numscalar("r(pb_model_failcode)", 0)
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
			real scalar normalized,
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
						B = delta * int_{cutoff}^{zbar} h0(z) dz / bw
						gamma = beta / (1 + delta)

					Use the mass equation to define:
						delta_U = B / (R * beta')
				*/

				R = intbasis(cutoff_est, zbar_est, K) / bw_est
				Rbeta = R * beta'

				if (Rbeta <= 0 | Rbeta >= . | B >= .) return

				delta = B / Rbeta

				if (1 + delta <= 0 | delta >= .) return

				q = (gamma - beta :/ (1 + delta))'

				/*
					delta = B / Rbeta
					ddelta/dbeta = -B * R / Rbeta^2
					ddelta/dB    = 1 / Rbeta
				*/
				ddelta_dbeta = -B * R / (Rbeta^2)

				ddelta_dtheta =
					ddelta_dbeta,
					J(1, Kb, 0),
					1/Rbeta

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
					normalized,
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
					normalized,
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
					normalized,
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
						normalized,
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
                else       G[7,.] = dx / (x * L)
            }
        }

        st_matrix(Gname, G)
    }

    return(out)
}



		end
