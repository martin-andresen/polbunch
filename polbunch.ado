		*! polbunch version date 20253001
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
			INITvals(string) ///
			notransform ///
			nopositiveshift ///
			nonormalize ///
			BOOTreps(integer 1) ///
			log ///
			constant ///
			nodots /// suppress dots for bootstrap progress
			notest /// do not test
			nosmallsample ///
			nobayes ///
			nozero ///
			saveunres(string) ///
			Bmodel ///
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
				
				if `bootreps'>0 {
					loc coeftabresults=c(coeftabresults)
					set coeftabresults off
				}
				if !inlist(`estimator',0,1,2,3) {
					noi di as error "Option estimator can take only values 0 (using data to the left only),  1 (no adjustment), 2 (Chetty et. al. adjustment) or 3 (theoretically consistent and efficient estimator)."
					exit 301
				}
				if "`test'" != "notest" {
					noi di as text "Note: model-restriction tests are disabled in the unified stacked branch."
					local test notest
				}
				tempvar touse
				marksample touse
				preserve
				drop if !`touse'
				
				if `bootreps'<0 {
					noi di as error "Option bootreps can only take values 0 (no inference, for outside bootstrap), 1 (analytic standard errors, the default) or a positive integer >1 (binned bootstrap)."
					exit 301
				}
				if `polynomial'<0 {
					noi di as error "Polynomial must be a positive integer"
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
				else { //collapse data
					loc z `varlist'
					tempvar bin y
					
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
					loc N=r(N)
					gen `bin'=ceil((`z'-`cutoff')/`bw')*`bw'+`cutoff'-`bw'/2
					
					collapse (count) `y'=`z', by(`bin')
					rename `bin' `z'
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
					replace `relbin' = -ceil((`cutoff' - `zright')/`bw') if `zright' <= `cutoff'
					replace `relbin' =  ceil((`zleft'  - `cutoff')/`bw') if `zleft'  >= `cutoff'

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
					replace `z' = (`z' - `cutoff') / `bw'
					local cutoff_est = 0
					local bw_est = 1
					local xscale = `bw'
				}
				else {
					local cutoff_est = `cutoff'
					local bw_est = `bw'
					local xscale = 1
				}
				
				tempname table
				mkmat `y' `z', matrix(`table')
				mat colnames `table'= freq `z'

				summarize `z', meanonly
				local zbar_est = r(max) + 0.5*`bw_est'
				
				tempvar side
				gen byte `side' = .
				replace `side' = -1 if `bunch' == 0 & `zright' <= `cutoff'
				replace `side' =  1 if `bunch' == 0 & `zleft'  >= `cutoff'

				count if `bunch' == 0 & missing(`side')
				if r(N) > 0 {
					noi di as error "Some non-excluded bins cannot be classified as left or right of cutoff."
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


				if `polynomial'>0 {
					loc nmiss=1
					while `nmiss'>0 {
						regress `y' 0.`dum'#(`rhsvars') 0.`dum' 1.`dum2'#(`rhsvars') 1.`dum2' if `bunch'==0, nocons
						loc nmiss=e(rank)<(`polynomial'+1)*2
						if `nmiss'>0 {
							loc note note
							loc polynomial=`polynomial'-1
							if `polynomial'<0 {
								noi di in red "Could not estimate separate polynomials on either side of the cutoff."
								exit 301
								}
							
							//NEW NAMES
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
						}
					if "`note'"=="note" {
						noi di as text "Note: Polynomial order lowered to `polynomial' because of multicollinearity problems with the specified polynomial."
					}
				} 
				else {
					loc unresmodel 0.`dum' 1.`dum' b0.`bunch'
				}
				
				//CHANGE NAMES for unrestricted model
				tempname bu Vu
				mat `bu'=e(b)
				
				loc unresnames `names' `names'
				loc unreseqnames `coleq0' `coleq1' 
				mat colnames `bu'=`unresnames'
				mat coleq `bu'=`unreseqnames'
				
				tempname h0coef h1coef
				matrix `h0coef' = J(1, `polynomial' + 1, .)
				matrix `h1coef' = J(1, `polynomial' + 1, .)

				forvalues j = 1/`=`polynomial'+1' {
					matrix `h0coef'[1,`j'] = `bu'[1,`j']
					matrix `h1coef'[1,`j'] = `bu'[1,`=`polynomial'+1+`j'']
				}

								
/*
				Simple delta initializer from one beta/gamma coefficient relation.
				h0coef/h1coef ordering:
					beta_1 ... beta_K beta_0
			*/

			scalar dstart = .

			if `estimator' == 2 {
				/*
					Prefer constant if usable, otherwise first polynomial coefficient.
				*/
				local j = `=`polynomial' + 1'

				if abs(`h1coef'[1,`j']) > 1e-12 {
					scalar dstart = `h0coef'[1,`j'] / `h1coef'[1,`j'] - 1
				}

				if missing(dstart) | dstart <= 0 {
					local j = 1
					if `polynomial' >= 1 & abs(`h1coef'[1,`j']) > 1e-12 {
						scalar dstart = `h0coef'[1,`j'] / `h1coef'[1,`j'] - 1
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
						if abs(`h0coef'[1,`j']) > 1e-12 & ///
						   `h1coef'[1,`j'] / `h0coef'[1,`j'] > 0 {
							scalar dstart = ///
								(`h1coef'[1,`j'] / `h0coef'[1,`j'])^(1/(`j' + 1)) - 1
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
				if `polynomial' >= 1 & abs(`h0coef'[1,1]) > 1e-12 {
					scalar dstart = exp( ///
						(`h1coef'[1,`=`polynomial' + 1'] - ///
						 `h0coef'[1,`=`polynomial' + 1']) / ///
						 `h0coef'[1,1] ///
					) - 1
				}
			}

			/*
				Fallback / bounds.
			*/
			if missing(dstart) | dstart <= 0 {
				scalar dstart = 0.05
			}

			if dstart > 2 {
				scalar dstart = 2
			}

			local dstart = scalar(dstart)
							
				//BOOTSTRAP SETUP
				if `bootreps'>1 {
					tempname p yorig
					if "`bayes'"=="nobayes" gen double `p'=`y'/`N'
					else {
						gen double `yorig'=`y'
						recast double `y'
					}
				}
						
				
				//ESTIMATION AND INFERENCE
				tempname b V bs tmpb bus Vus bb VV
				
				if inlist(`bootreps',0,1) loc stop=0
				else loc stop=`bootreps'
				forvalues s=0/`stop' {
					if `s'==1&"`dots'"!="nodots" nois _dots 0, title("Performing bootstrap repetitions...") reps(`bootreps')
					if `s'>0 { //resample outcome
						if "`bayes'"=="nobayes" {
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
						if `bootreps'==1 local vce vce
						bunch_profile `y' `z' `side' `bunch', ///
							estimator(`estimator') k(`polynomial') ///
							cutoff_orig(`cutoff_orig') bw_orig(`bw_orig') ///
							l(`L') h(`H') ///
							`log' `normalize' `vce' ///
							initdelta(`dstart') ///
							zbar_est(`zbar_est') ///
							zl_excl_orig(`zL_excl_orig') ///
							zh_excl_orig(`zH_excl_orig')
							
						matrix `b' = e(b)
						capture matrix `V' = e(V)
						
					////TRANSFORM ESTIMATES
					if "`transform'"!="notransform" {
						if `bootreps'!=1 loc nograd nograd
						quietly summarize `y' if `bunch' > 0, meanonly
						local Hstar_obs = r(sum)
						local taxopts
						if "`t0'" != "" & "`t1'" != "" {
							local taxopts t0(`t0') t1(`t1')
						}
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
							`log' ///
							`constant' ///
							`taxopts' ///
							`grad' ///
							`normalize' ///
							zbar(`zbar_est') ///
							massobs(`Hstar_obs') ///
							`bmodel'			
							
							matrix `b' = e(b)
							capture matrix `V' = e(V)
						}

						//BOOTSTRAP WRAPUP: COLLECT ESTIMATES
						if `bootreps'>1&`s'>0 mat `bs'=nullmat(`bs') \ e(b)
						if `s' > 0 noi _dots `s' 0
						
					}
					
						
					//bootstrap inference sunmmary & test
					if `bootreps'>1 {
						clear
						svmat `bs'
						corr _all, cov
						mat `V'=r(C)
						if `estimator'>0&"`test'"!="notest" {
							clear 
							svmat `bus'
							corr _all, cov
							mat `Vus'=r(C)
						}
					}
					
					//TEST RESTRICTIONS
					if `estimator'>0&"`test'"!="notest" {
						if `bootreps'>1 {
							ereturn post `bu' `Vus'		
							}
						else {
							est restore saveunres
							}
						testnl `teststr'
						loc chi2=r(chi2)
						loc p_mod=r(p)
						loc df=r(df)
						if "`saveunres'"!="" {
							ereturn local cmd="polbunch"
							est sto `saveunres'
						}
					}

				
				//POST RESULTS
				restore
				if `bootreps'>=1 eret post `b' `V', esample(`touse') depname(freq) obs(`N')
				else eret post `b', esample(`touse') obs(`N') depname(freq)
				if `estimator'>0&"`test'"!="notest"&`bootreps'>0 {
					estadd scalar chi2=`chi2'
					estadd scalar p_mod=`p_mod'
					estadd scalar df_mod=`df'
				}
				ereturn scalar polynomial=`polynomial'
				ereturn scalar bandwidth=`bw'
				ereturn scalar cutoff=`cutoff'
				ereturn scalar lower_limit=`zL_excl_orig'
				ereturn scalar upper_limit=`zH_excl_orig'
				ereturn local normalize="`normalize'"
				if `bootreps'>0 ereturn local cmd "polbunch"
				ereturn local cmdname "polbunch"
				ereturn local title 	"Polynomial bunching estimates"
				ereturn local cmdline 	"polbunch `0'"
				ereturn matrix table=`table'
				ereturn local binname "`z'"
				ereturn scalar bw=`bw'
				if `bootreps'>1 estadd local vcetype "bootstrap"
				if "`log'"=="log" ereturn scalar log=1
				else ereturn scalar log=0
	
				//Display results
				noi {
					di _newline
					di "`e(title)'"
					eret di
					if `estimator'>0&"`test'"!="notest"&`bootreps'>0 {
						di "Test of model assumptions: {col 42}Chi2(`df') test statistic {col 72}`: di %12.4f `chi2''"
						di "{col 42}p-value {col 72}`: di %12.4f `p_mod''"
						di "{hline 83}"
						}
					if inlist(`estimator',1,2) {
						di "Note: Estimator is not consistent with iso-elastic labor supply model and thus biased."
					}
					if "`constant'"!=""&"`t0'"!=""&"`t1'"!=""&"`transform'"!="notransform" {
						di "Note: Using the constant approximation to the density to calculate the elasticity may lead to bias."
					}
				}
			
			if `bootreps'>0 set coeftabresults `coeftabresults'
			}
				
		end
		
program define varcorrect, rclass
    syntax anything, [nosmallsample]

    quietly {
        gettoken y anything : anything

        tempvar res rss
        tempname g V
        local numbins = _N

        summarize `y', meanonly
        local N = r(sum)
        local k = e(df_m)

        predict `res', residuals
        return local xvars "`anything'"

        mata: st_matrix("`V'", ///
             varcorrect(st_data(., tokens(st_local("anything"))), ///
                          st_data(., "`y'"), ///
                          st_data(., "`res'"), 0))
        if "`smallsample'" != "nosmallsample" {
            matrix `V' = (`N'/(`N'-1)) * ///
                ((`N'*`numbins'-1)/(`N'*`numbins'-`k')) * `V'
        }

        return matrix V = `V'

        gen double `rss' = `y'*(`N'-`res')^2 + (`N'-`y')*(`res')^2
        summarize `rss', meanonly
        return local rss = r(sum)
    }
end

program define bunch_profile, eclass
    version 16.0

    syntax varlist(min=4 max=4 numeric) [if] [in] , ///
        ESTimator(integer) ///
        K(integer) ///
        CUTOFF_orig(real) ///
        BW_orig(real) ///
        L(integer) ///
        H(integer) ///
        zbar_est(real) ///
        zl_excl_orig(real) ///
        zh_excl_orig(real) ///
        [ nonormalize LOG VCE initdelta(real 0.05) ]

    gettoken yvar rest : varlist
    gettoken zvar rest : rest
    gettoken sidevar bunch : rest

    if !inlist(`estimator', 0, 1, 2, 3) {
        di as err "estimator() must be 0, 1, 2, or 3"
        exit 198
    }

   marksample touse, novarlist
   replace `touse' = 0 if missing(`yvar') | missing(`zvar') | missing(`bunch')

    local normalized0 = ("`normalize'" != "nonormalize")
    local islog0      = ("`log'"        != "")
    local dovar0      = ("`vce'"        != "")

    tempvar y_t z_t side_t bunch_t

    gen double `y_t'     = `yvar'   if `touse'
    gen double `z_t'     = `zvar'   if `touse'
    gen double `side_t'  = `sidevar' if `touse'
    gen double `bunch_t' = `bunch'  if `touse'

    tempname b V Gstack mustack minusmustack stackid

    mata: profile_run( ///
        "`y_t'", ///
        "`z_t'", ///
        "`side_t'", ///
        "`bunch_t'", ///
        `cutoff_orig', ///
        `bw_orig', ///
        `k', ///
        `estimator', ///
        `normalized0', ///
        `islog0', ///
        `zl_excl_orig', ///
        `zh_excl_orig', ///
        `zbar_est', ///
        `dovar0', ///
        `initdelta' ///
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

    if `dovar0' {
        ereturn local vcetype "Collapsed sandwich"

        matrix `Gstack'       = r_G_stack
        matrix `mustack'      = r_mu_stack
        matrix `minusmustack' = r_minus_mu_stack
        matrix `stackid'      = r_stack_id

        ereturn matrix G_stack        = `Gstack'
        ereturn matrix mu_stack       = `mustack'
        ereturn matrix minus_mu_stack = `minusmustack'
        ereturn matrix stack_id       = `stackid'
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
        LOW(integer) ///
        HIGH(integer) ///
        [ LOG CONSTANT T0(numlist min=1 max=1) T1(numlist min=1 max=1) NOGRAD nonormalize ZBAR(real 0) MASSOBS(real 0) BMODEL ]

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
    local hastax = "`t0'" != "" & "`t1'" != ""

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
    local hastax0     = ("`hastax'"   != "")
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
        `zhexcl' ///
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
            di as err "conformability error: colsof(G) != colsof(e(b))"
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
        ereturn local vcetype "Delta method"
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
    real colvector minus_mu
    real colvector stack_id
    real colvector fw_orig
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
    real scalar normalized,
    real scalar islog
)
{
    if (islog == 1) {
        return(ln(1 + delta))
    }

    if (normalized == 1) {
        return(delta * cutoff_orig / bw_orig)
    }

    return(delta * cutoff_orig)
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
    real scalar normalized,
    real scalar islog,
    real scalar dograd
)
{
    struct hcoef_out scalar out
    struct hcoef_out scalar hp, hm

    real scalar Kb, p, j
    real scalar scale, s, a
    real scalar eps, dp, dm
    real rowvector gp, gm

    Kb = K + 1

    out.gamma = J(1, Kb, .)

    if (dograd) {
        out.dgamma_dbeta  = J(Kb, Kb, 0)
        out.dgamma_ddelta = J(Kb, 1, 0)
    }
    else {
        out.dgamma_dbeta  = J(0, 0, .)
        out.dgamma_ddelta = J(0, 1, .)
    }

    if (estimator == 1) {
        out.gamma = beta

        if (dograd) {
            out.dgamma_dbeta = I(Kb)
        }
    }
    else if (estimator == 2) {
        out.gamma = beta :/ (1 + delta)

        if (dograd) {
            out.dgamma_dbeta  = I(Kb) / (1 + delta)
            out.dgamma_ddelta = -beta' / (1 + delta)^2
        }
    }
    else if (estimator == 3) {
        /*
            Generic affine map:
                h1(x) = scale * h0(a + s*x)

            beta ordering:
                beta[1]   coefficient on x^1
                ...
                beta[K]   coefficient on x^K
                beta[K+1] constant

            Level case:
                z0 = cutoff + (1+delta)*(z - cutoff)

                normalized x = (z-cutoff)/bw:
                    x0 = (1+delta)*x
                    scale = 1+delta
                    a = 0
                    s = 1+delta

                non-normalized z:
                    z0 = cutoff + (1+delta)*(z-cutoff)
                       = (1-s)*cutoff + s*z
                    scale = 1+delta
                    a = (1-s)*cutoff
                    s = 1+delta

            Log/proportional case:
                z0 = (1+delta)*z

                normalized x = (z-cutoff)/bw:
                    x0 = ((1+delta)*(cutoff + bw*x) - cutoff)/bw
                       = delta*cutoff/bw + (1+delta)*x
                    scale = 1+delta
                    a = delta*cutoff/bw
                    s = 1+delta

                non-normalized z:
                    z0 = (1+delta)*z
                    scale = 1+delta
                    a = 0
                    s = 1+delta
        */

        scale = 1 + delta
        s     = 1 + delta

        if (islog == 0) {
            if (normalized == 1) {
                a = 0
            }
            else {
                a = (1 - s) * cutoff_orig
            }
        }
        else {
            if (normalized == 1) {
                a = delta * cutoff_orig / bw_orig
            }
            else {
                a = 0
            }
        }

        out.gamma = J(1, Kb, 0)

        if (dograd) {
            out.dgamma_dbeta = J(Kb, Kb, 0)
        }

        /*
            For p = 1,...,K:
                coeff on x^p in scale * beta_j * (a+s*x)^j
                equals scale * beta_j * comb(j,p) * a^(j-p) * s^p
        */
        for (p = 1; p <= K; p++) {
            for (j = p; j <= K; j++) {
                out.gamma[p] = out.gamma[p] +
                    beta[j] * scale * comb(j,p) * a^(j-p) * s^p

                if (dograd) {
                    out.dgamma_dbeta[p,j] =
                        scale * comb(j,p) * a^(j-p) * s^p
                }
            }
        }

        /*
            Constant:
                scale * beta_0 + scale * sum_j beta_j * a^j
        */
        out.gamma[Kb] = scale * beta[Kb]

        if (dograd) {
            out.dgamma_dbeta[Kb,Kb] = scale
        }

        for (j = 1; j <= K; j++) {
            out.gamma[Kb] = out.gamma[Kb] + scale * beta[j] * a^j

            if (dograd) {
                out.dgamma_dbeta[Kb,j] = scale * a^j
            }
        }

        if (dograd) {
            /*
                Finite-difference derivative wrt delta. This keeps level/log
                and normalized/non-normalized cases consistent without deriving
                four separate analytic derivatives.
            */
            eps = max((1e-6, abs(delta)*1e-5))
            dp  = delta + eps
            dm  = max((delta - eps, 1e-10))

            hp = h1coef_map(beta, dp, estimator, K,
                cutoff_orig, bw_orig, normalized, islog, 0)

            hm = h1coef_map(beta, dm, estimator, K,
                cutoff_orig, bw_orig, normalized, islog, 0)

            gp = hp.gamma
            gm = hm.gamma

            out.dgamma_ddelta = ((gp - gm) / (dp - dm))'
        }
    }
    else {
        _error(3498, "h1coef_map only handles estimators 1, 2, and 3")
    }

    return(out)
}

       
// design row transformation for h1, consistent with h1coef_map()
struct hdesign_out scalar h1design23(
    real scalar delta,
    real colvector zR,
    real scalar cutoff_orig,
    real scalar bw_orig,
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar dograd
)
{
    struct hdesign_out scalar out
    struct hcoef_out scalar h1, hp, hm

    real scalar eps, dp, dm
    real matrix Xbase, Xp, Xm

    Xbase = pbasis(zR, K)

    h1 = h1coef_map(
        J(1, K+1, 0),
        delta,
        estimator,
        K,
        cutoff_orig,
        bw_orig,
        normalized,
        islog,
        1
    )

    /*
        If gamma = A(delta) * beta', then fitted h1 rows are:
            Xbase * A(delta) * beta'
    */
    out.X = Xbase * h1.dgamma_dbeta

    if (dograd) {
        eps = max((1e-6, abs(delta)*1e-5))
        dp  = delta + eps
        dm  = max((delta - eps, 1e-10))

        hp = h1coef_map(
            J(1, K+1, 0),
            dp,
            estimator,
            K,
            cutoff_orig,
            bw_orig,
            normalized,
            islog,
            1
        )

        hm = h1coef_map(
            J(1, K+1, 0),
            dm,
            estimator,
            K,
            cutoff_orig,
            bw_orig,
            normalized,
            islog,
            1
        )

        Xp = Xbase * hp.dgamma_dbeta
        Xm = Xbase * hm.dgamma_dbeta

        out.dXddelta = (Xp - Xm) / (dp - dm)
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
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar zbar_est
)
{
    real scalar cutoff_est, bw_est, r
    real rowvector R

    if (normalized == 1) {
        cutoff_est = 0
        bw_est     = 1
    }
    else {
        cutoff_est = cutoff_orig
        bw_est     = bw_orig
    }

    if (estimator == 2) {
        /* Chetty restriction: bw * B = delta * int_{zstar}^{zbar} h0(z) dz */
        R = (delta / bw_est) * intbasis(cutoff_est, zbar_est, K)
    }
    else if (estimator == 3) {
        /* Theoretically consistent restriction: bw * B = int_{zstar}^{zstar+r(delta)} h0(z) dz */
        r = response_length(delta, cutoff_orig, bw_orig, normalized, islog)
        R = intbasis(cutoff_est, cutoff_est + r, K) / bw_est
    }
    else {
        _error(3498, "bmodel_row23 only handles estimators 2 and 3")
    }

    return(R)
}

real rowvector d_bmodel_row23_ddelta(
    real scalar delta,
    real scalar cutoff_orig,
    real scalar bw_orig,
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar zbar_est
)
{
    real scalar eps, dp, dm
    real rowvector Rp, Rm

    eps = max((1e-6, abs(delta)*1e-5))
    dp  = delta + eps
    dm  = max((delta - eps, 1e-10))

    Rp = bmodel_row23(dp, cutoff_orig, bw_orig, K, estimator, normalized, islog, zbar_est)
    Rm = bmodel_row23(dm, cutoff_orig, bw_orig, K, estimator, normalized, islog, zbar_est)

    return((Rp - Rm)/(dp - dm))
}

// cf_mass_row23() returns R(delta) such that cf_mass(beta,delta) = R(delta) * beta'.
// cf_mass = integral left part under h0 + right part under h1.
real rowvector cf_mass_row23(
    real scalar delta,
    real colvector z,
    real scalar cutoff_orig,
    real scalar bw_orig,
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar zL_excl_orig,
    real scalar zH_excl_orig
)
{
    real scalar ex_lo, ex_hi, cutoff_est, bw_est
    real rowvector Rlo, Rhi
    struct hcoef_out scalar h1map

	if (normalized == 1) {
		cutoff_est = 0
		bw_est     = 1
		ex_lo      = (zL_excl_orig - cutoff_orig) / bw_orig
		ex_hi      = (zH_excl_orig - cutoff_orig) / bw_orig
	}
	else {
		cutoff_est = cutoff_orig
		bw_est     = bw_orig
		ex_lo      = zL_excl_orig
		ex_hi      = zH_excl_orig
	}
	
    Rlo = intbasis(ex_lo, cutoff_est, K)

    if (estimator == 2) {
        Rhi = intbasis(cutoff_est, ex_hi, K) / (1 + delta)
    }
    else if (estimator == 3) {
        h1map = h1coef_map(
			J(1, K+1, 0),
			delta,
			estimator,
			K,
			cutoff_orig,
			bw_orig,
			normalized,
			islog,
			1
		)
        Rhi = intbasis(cutoff_est, ex_hi, K) * h1map.dgamma_dbeta
    }
    else {
        _error(3498, "cf_mass_row23 only handles estimators 2 and 3")
    }

    return((Rlo + Rhi) / bw_est)
}

real rowvector d_cf_mass_row23_ddelta(
    real scalar delta,
    real colvector z,
    real scalar cutoff_orig,
    real scalar bw_orig,
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar zL_excl_orig,
    real scalar zH_excl_orig
)
{
    real scalar eps, dp, dm
    real rowvector Rp, Rm

    eps = max((1e-6, abs(delta)*1e-5))
    dp  = delta + eps
    dm  = max((delta - eps, 1e-10))

    Rp = cf_mass_row23(dp, z, cutoff_orig, bw_orig, K, estimator, normalized, islog, zL_excl_orig, zH_excl_orig)
    Rm = cf_mass_row23(dm, z, cutoff_orig, bw_orig, K, estimator, normalized, islog, zL_excl_orig, zH_excl_orig)

    return((Rp - Rm)/(dp - dm))
}

real rowvector h0_excluded_row(
    real scalar cutoff_orig,
    real scalar bw_orig,
    real scalar K,
    real scalar normalized,
    real scalar zL_excl_orig,
    real scalar zH_excl_orig
)
{
    real scalar ex_lo, ex_hi, bw_est

    if (normalized == 1) {
        bw_est = 1
        ex_lo  = (zL_excl_orig - cutoff_orig) / bw_orig
        ex_hi  = (zH_excl_orig - cutoff_orig) / bw_orig
    }
    else {
        bw_est = bw_orig
        ex_lo  = zL_excl_orig
        ex_hi  = zH_excl_orig
    }

    return(intbasis(ex_lo, ex_hi, K) / bw_est)
}
// -----------------------------------------------------------------------------
// Unified stacked design/profile objective for estimators 0/1/2/3
// -----------------------------------------------------------------------------

real rowvector cf_mass_row(
    real scalar delta,
    real scalar cutoff_orig,
    real scalar bw_orig,
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar zL_excl_orig,
    real scalar zH_excl_orig,
    real scalar ntheta
)
{
    real scalar ex_lo, ex_hi, cutoff_est, bw_est, Kb
    real rowvector R, Rlo, Rhi
    struct hcoef_out scalar h1map

    Kb = K + 1

    if (normalized == 1) {
        cutoff_est = 0
        bw_est     = 1
        ex_lo      = (zL_excl_orig - cutoff_orig) / bw_orig
        ex_hi      = (zH_excl_orig - cutoff_orig) / bw_orig
    }
    else {
        cutoff_est = cutoff_orig
        bw_est     = bw_orig
        ex_lo      = zL_excl_orig
        ex_hi      = zH_excl_orig
    }

    R = J(1, ntheta, 0)

    if (estimator == 0) {
        Rlo = intbasis(ex_lo, cutoff_est, K) / bw_est
        Rhi = intbasis(cutoff_est, ex_hi, K) / bw_est
        R[1, 1..Kb] = Rlo
        R[1, (Kb+1)..(2*Kb)] = Rhi
    }
    else if (estimator == 1) {
        R[1, 1..Kb] = intbasis(ex_lo, ex_hi, K) / bw_est
    }
    else if (estimator == 2) {
        Rlo = intbasis(ex_lo, cutoff_est, K) / bw_est
        Rhi = intbasis(cutoff_est, ex_hi, K) / ((1 + delta) * bw_est)
        R[1, 1..Kb] = Rlo + Rhi
    }
    else if (estimator == 3) {
        Rlo = intbasis(ex_lo, cutoff_est, K) / bw_est
        h1map = h1coef_map(
            J(1, K+1, 0),
            delta,
            estimator,
            K,
            cutoff_orig,
            bw_orig,
            normalized,
            islog,
            1
        )
        Rhi = (intbasis(cutoff_est, ex_hi, K) * h1map.dgamma_dbeta) / bw_est
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
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar zL_excl_orig,
    real scalar zH_excl_orig,
    real scalar ntheta
)
{
    real scalar eps, dp, dm
    real rowvector Rp, Rm

    if (estimator == 0 | estimator == 1) {
        return(J(1, ntheta, 0))
    }

    eps = max((1e-6, abs(delta)*1e-5))
    dp  = delta + eps
    dm  = max((delta - eps, 1e-10))

    Rp = cf_mass_row(dp, cutoff_orig, bw_orig, K, estimator, normalized, islog, zL_excl_orig, zH_excl_orig, ntheta)
    Rm = cf_mass_row(dm, cutoff_orig, bw_orig, K, estimator, normalized, islog, zL_excl_orig, zH_excl_orig, ntheta)

    return((Rp - Rm)/(dp - dm))
}

real rowvector bmodel_row(
    real scalar delta,
    real scalar cutoff_orig,
    real scalar bw_orig,
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar zbar_est,
    real scalar ntheta
)
{
    real scalar cutoff_est, bw_est, r
    real rowvector R

    R = J(1, ntheta, 0)

    if (estimator == 0 | estimator == 1) {
        return(R)
    }

    if (normalized == 1) {
        cutoff_est = 0
        bw_est     = 1
    }
    else {
        cutoff_est = cutoff_orig
        bw_est     = bw_orig
    }

    if (estimator == 2) {
        R[1, 1..(K+1)] = (delta / bw_est) * intbasis(cutoff_est, zbar_est, K)
    }
    else if (estimator == 3) {
        r = response_length(delta, cutoff_orig, bw_orig, normalized, islog)
        R[1, 1..(K+1)] = intbasis(cutoff_est, cutoff_est + r, K) / bw_est
    }
    else {
        _error(3498, "bmodel_row only handles estimators 0, 1, 2, and 3")
    }

    return(R)
}

real rowvector d_bmodel_row_ddelta(
    real scalar delta,
    real scalar cutoff_orig,
    real scalar bw_orig,
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar zbar_est,
    real scalar ntheta
)
{
    real scalar eps, dp, dm
    real rowvector Rp, Rm

    if (estimator == 0 | estimator == 1) {
        return(J(1, ntheta, 0))
    }

    eps = max((1e-6, abs(delta)*1e-5))
    dp  = delta + eps
    dm  = max((delta - eps, 1e-10))

    Rp = bmodel_row(dp, cutoff_orig, bw_orig, K, estimator, normalized, islog, zbar_est, ntheta)
    Rm = bmodel_row(dm, cutoff_orig, bw_orig, K, estimator, normalized, islog, zbar_est, ntheta)

    return((Rp - Rm)/(dp - dm))
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
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar zL_excl_orig,
    real scalar zH_excl_orig,
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

        Xcf = cf_mass_row(delta, cutoff_orig, bw_orig, K, estimator, normalized, islog, zL_excl_orig, zH_excl_orig, ntheta)
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

        Xcf = cf_mass_row(delta, cutoff_orig, bw_orig, K, estimator, normalized, islog, zL_excl_orig, zH_excl_orig, ntheta)
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
        h1 = h1design23(delta, zR, cutoff_orig, bw_orig, K, estimator, normalized, islog, dograd)

        ntheta = Kb
        Xcf    = cf_mass_row(delta, cutoff_orig, bw_orig, K, estimator, normalized, islog, zL_excl_orig, zH_excl_orig, ntheta)
        Xbmod  = bmodel_row(delta, cutoff_orig, bw_orig, K, estimator, normalized, islog, zbar_est, ntheta)
        Xmass  = Xcf + Xbmod

        out.X = XL \ h1.X \ Xmass

        if (dograd) {
            dX0 = J(rows(XL), Kb, 0)
            dXcf   = d_cf_mass_row_ddelta(delta, cutoff_orig, bw_orig, K, estimator, normalized, islog, zL_excl_orig, zH_excl_orig, ntheta)
            dXbmod = d_bmodel_row_ddelta(delta, cutoff_orig, bw_orig, K, estimator, normalized, islog, zbar_est, ntheta)
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
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar zL_excl_orig,
    real scalar zH_excl_orig,
    real scalar zbar_est
)
{
    real colvector ystack, theta
    struct design_out scalar D

    ystack = make_ystack(y, side, bunch, estimator, Hstar_obs)

    D = make_design(delta, z, side, bunch, cutoff_orig, bw_orig, K, estimator, normalized, islog, zL_excl_orig, zH_excl_orig, zbar_est, 0)

    theta = qrsolve(D.X, ystack)
    return(theta')
}

real scalar profQ(
    real scalar lndelta,
    real colvector y,
    real colvector z,
    real colvector side,
    real colvector bunch,
    real scalar Hstar_obs,
    real scalar cutoff_orig,
    real scalar bw_orig,
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar zL_excl_orig,
    real scalar zH_excl_orig,
    real scalar zbar_est
)
{
    real scalar delta
    real colvector ystack, theta, resid
    struct design_out scalar D

    if (estimator == 0 | estimator == 1) {
        return(0)
    }

    delta = exp(lndelta)
    ystack = make_ystack(y, side, bunch, estimator, Hstar_obs)

    D = make_design(delta, z, side, bunch, cutoff_orig, bw_orig, K, estimator, normalized, islog, zL_excl_orig, zH_excl_orig, zbar_est, 0)

    theta = qrsolve(D.X, ystack)
    resid = ystack - D.X * theta

    return(quadcross(resid, resid))
}

struct stack23_out scalar profile_stack(
    real colvector y,
    real colvector z,
    real colvector side,
    real colvector bunch,
    real rowvector theta,
    real scalar delta,
    real scalar Hstar_obs,
    real scalar cutoff_orig,
    real scalar bw_orig,
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar zL_excl_orig,
    real scalar zH_excl_orig,
    real scalar zbar_est,
    real scalar dograd
)
{
    struct stack23_out scalar out
    struct design_out scalar D

    real scalar i, idx, nleft, nright, nout, Kb
    real rowvector beta

    Kb = K + 1

    out.ystack = make_ystack(y, side, bunch, estimator, Hstar_obs)

    D = make_design(delta, z, side, bunch, cutoff_orig, bw_orig, K, estimator, normalized, islog, zL_excl_orig, zH_excl_orig, zbar_est, dograd)

    out.X = D.X

    if (estimator == 0 | estimator == 1) {
        out.mu = D.X * theta'
        out.G  = D.X
    }
    else {
        beta   = theta[1..Kb]
        out.mu = D.X * beta'

        if (dograd) {
            out.G = D.X, D.dXddelta * beta'
        }
        else {
            out.G = J(0, 0, .)
        }
    }

    out.minus_mu = out.ystack - out.mu
    out.stack_id = J(rows(y), 1, .)

    if (estimator == 1) {
        nout = sum(bunch :== 0)
        idx = 0
        for (i = 1; i <= rows(y); i++) {
            if (bunch[i] == 0) {
                idx = idx + 1
                out.stack_id[i] = idx
            }
        }
        for (i = 1; i <= rows(y); i++) {
            if (bunch[i] > 0) {
                out.stack_id[i] = nout + 1
            }
        }
    }
    else {
        nleft  = sum((bunch :== 0) :& (side :== -1))
        nright = sum((bunch :== 0) :& (side :==  1))

        idx = 0
        for (i = 1; i <= rows(y); i++) {
            if (bunch[i] == 0 & side[i] == -1) {
                idx = idx + 1
                out.stack_id[i] = idx
            }
        }
        for (i = 1; i <= rows(y); i++) {
            if (bunch[i] == 0 & side[i] == 1) {
                idx = idx + 1
                out.stack_id[i] = idx
            }
        }
        for (i = 1; i <= rows(y); i++) {
            if (bunch[i] > 0) {
                out.stack_id[i] = nleft + nright + 1
            }
        }
    }

    out.fw_orig = y
    _editmissing(out.fw_orig, 0)

    return(out)
}

// -----------------------------------------------------------------------------
// Variance correction
// -----------------------------------------------------------------------------

real matrix varcorrect(real matrix X, real matrix fw, real matrix e, real scalar addcons)
{
    real scalar B, N, i
    real matrix meat, bread

    B = rows(fw)
    N = sum(fw)

    if (addcons == 1) {
        X = X, J(rows(X), 1, 1)
    }

    meat = J(cols(X), cols(X), 0)

    for (i = 1; i <= B; i++) {
        e[i] = e[i] + N
        meat = meat + fw[i] * (X' * e * e' * X)
        e[i] = e[i] - N
    }

    bread = invsym(quadcross(X, X) :* N)

    return(bread * meat * bread)
}

real matrix varcorrect_collapsed(
    real matrix G_stack,
    real colvector fw_orig,
    real colvector stack_id,
    real colvector e_stack,
    real scalar addcons
)
{
    real scalar B_orig, M, N, i, s
    real colvector e_i
    real matrix G, meat, bread

    G = G_stack

    if (addcons == 1) {
        G = G, J(rows(G), 1, 1)
    }

    B_orig = rows(fw_orig)
    M      = rows(G)
    N      = sum(fw_orig)

    meat = J(cols(G), cols(G), 0)

    for (i = 1; i <= B_orig; i++) {
        s = stack_id[i]

        if (s < .) {
            e_i = e_stack
            e_i[s] = e_i[s] + N

            meat = meat + fw_orig[i] * (G' * e_i * e_i' * G)
        }
    }

    bread = pinv(quadcross(G, G) :* N)
    return(bread * meat * bread)
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

/*
real scalar eresp(real scalar B, real scalar tau, real matrix cf, real scalar bw, real scalar xscale)
{
    real matrix cfpoly, integral, roots, realroots, out

    // Input cf is [b1, b2, ..., bK, b0]. Mata poly* functions need [b0, b1, ..., bK].
    if (cols(cf) == 1) {
        cfpoly = cf
    }
    else {
        cfpoly = cf[cols(cf)], cf[1..cols(cf)-1]
    }

    if (cols(cfpoly) == 1) {
        if (cfpoly[1] <= 0) return(.)
        return(((B * bw) / cfpoly[1]) * xscale)
    }

    integral = polyinteg(cfpoly, 1)
    integral[1] = -polyeval(integral, tau) - B*bw

    roots = polyroots(integral)
    realroots = Re(select(roots, Im(roots) :== 0))
    out = sort(select(realroots, realroots :> tau)', 1)'

    if (cols(out) == 0) return(.)
    else return((out[1] - tau) * xscale)
}
*/


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
	real scalar zH_excl_orig
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
            RB = bmodel_row23(delta, cutoff_orig, bw_orig, K, estimator, normalized, islog, zbar_est)

            B = RB * beta'
            b[1,oB] = B

            if (dograd) {
                G[oB, ibeta] = RB
                dRB = d_bmodel_row23_ddelta(delta, cutoff_orig, bw_orig, K, estimator, normalized, islog, zbar_est)
                G[oB, idelta] = dRB * beta'
            }
        }
        else if (Btype2 == 1) {
            /* Estimator 2, Chetty-style reduced-form B: B_reduced = Hstar_obs - int_{lower}^{upper} h0(z) dz / bw */
            RB = h0_excluded_row(cutoff_orig, bw_orig, K, normalized, zL_excl_orig, zH_excl_orig)

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
        RB = bmodel_row23(delta, cutoff_orig, bw_orig, K, estimator, normalized, islog, zbar_est)

        B = RB * beta'
        b[1,oB] = B

        if (dograd) {
            G[oB, ibeta] = RB
            dRB = d_bmodel_row23_ddelta(delta, cutoff_orig, bw_orig, K, estimator, normalized, islog, zbar_est)
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
    real scalar Hstar_obs, cutoff_orig, bw_orig, K
    real scalar estimator, normalized, islog
    real scalar zL_excl_orig, zH_excl_orig, zbar_est

    Hstar_obs    = pars[1]
    cutoff_orig  = pars[2]
    bw_orig      = pars[3]
    K            = pars[4]
    estimator    = pars[5]
    normalized   = pars[6]
    islog        = pars[7]
    zL_excl_orig = pars[8]
    zH_excl_orig = pars[9]
    zbar_est     = pars[10]

    val = profQ(
        p[1],
        y,
        z,
        side,
        bunch,
        Hstar_obs,
        cutoff_orig,
        bw_orig,
        K,
        estimator,
        normalized,
        islog,
        zL_excl_orig,
        zH_excl_orig,
        zbar_est
    )
}

void profile_run(
    string scalar yvar,
    string scalar zvar,
    string scalar sidevar,
    string scalar bunchvar,
    real scalar cutoff_orig,
    real scalar bw_orig,
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar zL_excl_orig,
    real scalar zH_excl_orig,
    real scalar zbar_est,
    real scalar dovar,
    real scalar initdelta
)
{
    real colvector y, z, side, bunch
    real scalar Kb, Hstar_obs, lndelta_hat, delta_hat, gi
    real rowvector theta_hat, beta_hat, b, pars, phat, dgrid, qgrid
    real matrix Vout
    transmorphic S
    struct stack23_out scalar st

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
            K,
            estimator,
            normalized,
            islog,
            zL_excl_orig,
            zH_excl_orig,
            zbar_est
        )
        b = theta_hat
    }
    else if (estimator == 2 | estimator == 3) {
        pars = (
            Hstar_obs,
            cutoff_orig,
            bw_orig,
            K,
            estimator,
            normalized,
            islog,
            zL_excl_orig,
            zH_excl_orig,
            zbar_est
        )

        S = optimize_init()
        optimize_init_evaluator(S, &profQ_opt())
        optimize_init_evaluatortype(S, "d0")
        optimize_init_which(S, "min")
        optimize_init_params(S, ln(initdelta))
        optimize_init_conv_maxiter(S, 200)

        optimize_init_argument(S, 1, y)
        optimize_init_argument(S, 2, z)
        optimize_init_argument(S, 3, side)
        optimize_init_argument(S, 4, bunch)
        optimize_init_argument(S, 5, pars)

        dgrid = (0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50)
        qgrid = J(1, cols(dgrid), .)

        for (gi = 1; gi <= cols(dgrid); gi++) {
            qgrid[gi] = profQ(
                ln(dgrid[gi]),
                y,
                z,
                side,
                bunch,
                Hstar_obs,
                cutoff_orig,
                bw_orig,
                K,
                estimator,
                normalized,
                islog,
                zL_excl_orig,
                zH_excl_orig,
                zbar_est
            )
        }

        st_matrix("debug_delta_grid", dgrid)
        st_matrix("debug_Q_grid", qgrid)

        phat         = optimize(S)
        lndelta_hat = phat[1]
        delta_hat   = exp(lndelta_hat)

        beta_hat = profTheta(
            delta_hat,
            y,
            z,
            side,
            bunch,
            Hstar_obs,
            cutoff_orig,
            bw_orig,
            K,
            estimator,
            normalized,
            islog,
            zL_excl_orig,
            zH_excl_orig,
            zbar_est
        )

        theta_hat = beta_hat, delta_hat
        b = theta_hat
    }
    else {
        _error(3498, "profile_run only handles estimators 0, 1, 2, and 3")
    }

    if (dovar == 1) {
        st = profile_stack(
            y,
            z,
            side,
            bunch,
            theta_hat,
            delta_hat,
            Hstar_obs,
            cutoff_orig,
            bw_orig,
            K,
            estimator,
            normalized,
            islog,
            zL_excl_orig,
            zH_excl_orig,
            zbar_est,
            1
        )
        Vout = varcorrect_collapsed(st.G, st.fw_orig, st.stack_id, st.minus_mu, 0)
    }
    else {
        Vout = J(cols(b), cols(b), .)
    }

    st_matrix("r_b_profile", b)
    st_matrix("r_V_profile", Vout)

    if (dovar == 1) {
        st_matrix("r_G_stack", st.G)
        st_matrix("r_mu_stack", st.mu)
        st_matrix("r_minus_mu_stack", st.minus_mu)
        st_matrix("r_stack_id", st.stack_id)
    }
}


end
