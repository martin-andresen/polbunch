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
				if !inlist(`estimator',0,1,2,3,4) {
					noi di as error "Option estimator can take only values 0 (using data to the left only),  1 (no adjustment), 2 (Chetty et. al. adjustment) or 3 (theoretically consistent and efficient estimator) or 4 (Saez trapezoid approximation)."
					exit 301
				}
				tempvar touse
				marksample touse
				preserve
				drop if !`touse'
				
				if `bootreps'<0 {
					noi di as error "Option bootreps can only take values 0 (no inference), 1 (analytic standard errors, the default) or a positive integer >0 (binned bootstrap)."
					exit 301
				}
				if `polynomial'<0 {
					noi di as error "Polynomial must be a positive integer"
					exit 301
				}
				
				
				// check varlist vs bw opts
				loc nvars: word count `varlist' 
				if (`nvars'==1&"`bw'"=="")|(`nvars'==2&"`bw'"!="") {
					noi di as error "Varlist must either contain 1 variable (earnings z) and option bw be specified (individual level data) or 2 variables (frequency y and earnings bin z) and option bandwith not be specified (pre-binned data)."
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
				

			// limits option
			if "`limits'"=="" {
				loc L = 1
				loc H = 0
			}
			else gettoken L H : limits

			// For pre-binned data, limits are numbers of whole bins below/above cutoff.
			// Do not split bins.  relbin = -1 is closest bin below cutoff;
			// relbin = 0 is first bin at or above cutoff.
			tempvar zleft zright relbin

			if `nvars' == 2 {
				gen double `zleft'  = `z' - `bw'/2
				gen double `zright' = `z' + `bw'/2

				gen int `relbin' = floor((`zleft' - `cutoff') / `bw') + 1

				summarize `z' if inrange(`relbin', -`L' + 1, `H'), meanonly
				local zL = r(min) - `bw'/2
				local zH = r(max) + `bw'/2
			}
			else {
				local zH = `cutoff' + `H'*`bw'
				local zL = `cutoff' - `L'*`bw'
			}
							
				//Put bins in e(table)
				tempname table
				mkmat `y' `z', matrix(`table')
				mat colnames `table'= freq `z'
				
				//check if there are people on either side
				count if `z'<`zL'
				if r(N)==0 {
					noi di as error "No individuals in sample allocates below the excluded region."
					exit 301
					}
				
				count if `z'>`zH'
				if r(N)==0 {
					if (`estimator'>0) {
						noi di as error "No individuals in sample allocates above the excluded region."
						exit 301
					}
					}
				
				//Reset pol for saez estimator
				if `estimator'==4 loc polynomial=0
				
				//NORMALIZE Z
				tempvar z_orig
				gen double `z_orig'=`z' 
				loc cutoff_orig = `cutoff'
				loc bw_orig = `bw'
				
				if "`normalize'" != "nonormalize" {
					replace `z' = (`z' - `cutoff') / `bw'
					local cutoff = 0
					local bw = 1
				}
				else {
					local cutoff= `cutoff'
					local bw = `bw'
				}
				
				//NAMES
				if `polynomial'>0 {
					forvalues i=1/`polynomial' { 
						if "`rhsvars'"=="" loc rhsvars  c.`z'
						else loc rhsvars `rhsvars'##c.`z'
						loc coleq0 `coleq0' h0
						loc coleq1 `coleq1' h1
						}
				
					loc coleq0 `coleq0' h0
					loc coleq1 `coleq1' h1
					
					fvexpand `rhsvars'
					loc names `r(varlist)' _cons
				}
				else {
					loc rhsvars 
					loc names _cons
					
					loc coleq0 h0
					loc coleq1 h1
				}
				
				//gen dummies
				tempvar fw dupe dum dum2 bunch cons
				gen byte `dum'=`z'>`cutoff'
				gen byte `dum2'=`dum'
				if `nvars' == 2 {
					egen `bunch' = group(`z') if inrange(`relbin', -`L' + 1, `H')
				}
				else {
					egen `bunch' = group(`z') if inrange(`z', `zL', `zH')
				}
				replace `bunch' = 0 if missing(`bunch')
				gen `cons'=1
				count
				loc numbins=r(N)
				

				//Evaluate multicollinearity
				if `polynomial'>0 {
					loc nmiss=1
					while `nmiss'>0 {
						regress `y' 0.`dum'#(`rhsvars') 0.`dum' 1.`dum2'#(`rhsvars') 1.`dum2' b0.`bunch', nocons
						loc nmiss=e(rank)<(`polynomial'+1)*2+`H'+`L'
						if `nmiss'>0 {
							loc note note
							loc polynomial=`polynomial'-1
							if `polynomial'<0 {
								noi di in red "Could not estimate separate polynomials on either side of the cutoff."
								exit 301
								}
							
							//NEW NAMES
							loc coleq0 
							loc coleq1 
							forvalues i=1/`polynomial' { 
								if `i'==1 loc rhsvars  c.`z'
								else loc rhsvars `rhsvars'##c.`z'
								loc coleq0 `coleq0' h0
								loc coleq1 `coleq1' h1
							}
					
							loc coleq0 `coleq0' h0
							loc coleq1 `coleq1' h1
							
							fvexpand `rhsvars'
							loc names `r(varlist)' _cons
							}
						}
					if "`note'"=="note" {
						noi di as text "Note: Polynomial order lowered to `polynomial' because of multicollinearity problems with the specified polynomial."
					}
				}
					
				//NAMES and model string FOR UNRESTRICTED MODEL (as benchmark or main model if estimator==0)
				forvalues bval=1/`=`H'+`L'' {
					loc bunchvars `bunchvars' `bval'.`bunch'
					}
				if `polynomial' == 0 {
					local unresmodel `y' 0.`dum' 1.`dum2' `bunchvars'
				}
				else {
					local unresmodel `y' 0.`dum'#(`rhsvars') 0.`dum' ///
						1.`dum2'#(`rhsvars') 1.`dum2' `bunchvars'
				}

				foreach l in b g {
					forvalues k = 1/`polynomial' {
						local newnames `newnames' /`l'`k'
						if `estimator' == 1 & "`l'" == "b" {
							local namesest1 `namesest1' /b`k'
						}
					}

					local newnames `newnames' /`l'0
					if `estimator' == 1 & "`l'" == "b" {
						local namesest1 `namesest1' /b0
					}
				}
				
				forvalues bval=1/`=`H'+`L'' {
					loc newnames `newnames' /bunch`bval'
					loc namesest1 `namesest1' /bunch`bval'
					}

				if `bootreps'>1 {
					tempname p yorig
					if "`bayes'"=="nobayes" gen double `p'=`y'/`N'
					else {
						gen double `yorig'=`y'
						recast double `y'
					}
				}

				// NAMES for final model
				if "`transform'" != "notransform" {
					if `estimator' != 2 {
						local names `names' `names' number_bunchers excess_mass shift marginal_response
					}
					else {
						local names `names' `names' delta number_bunchers excess_mass shift marginal_response
					}

					if "`t0'" != "" & "`t1'" != "" {
						local names `names' elasticity
					}
				}
				else {
					tokenize `names'
					local names
					local bnames

					forvalues i = 1/`L' {
						local bnames `bnames' number_bunchers:`=`L'-`i'+1'.below
					}

					if `H' > 0 {
						forvalues i = 1/`H' {
							local bnames `bnames' number_bunchers:`i'.above
						}
					}

					if `estimator' > 1 {
						if `estimator' == 2 {
							if "`positiveshift'" != "nopositiveshift" local names `names' delta:lndelta
							else local names `names' delta:delta
						}
						else {
							if "`positiveshift'" != "nopositiveshift" local names `names' shift:lnshift
							else local names `names' shift:shift
						}

						if `polynomial' > 0 {
							forvalues k = 1/`polynomial' {
								local names `names' h0:``k''
							}
						}
						local names `names' h0:_cons
						local names `bnames' `names'
					}
					else {
						// estimator 0 and 1: h0 coefficients
						if `polynomial' > 0 {
							forvalues k = 1/`polynomial' {
								local names `names' h0:``k''
							}
						}
						local names `names' h0:_cons

						// estimator 0 only: h1 coefficients
						if `estimator' == 0 {
							if `polynomial' > 0 {
								forvalues k = 1/`polynomial' {
									local names `names' h1:``k''
								}
							}
							local names `names' h1:_cons
						}

						local names `names' `bnames'
					}
					}							
							
				//BUILD TEST RESTRICITONS
				if "`test'"!="notest" {

				// Build symbolic restrictions g_k = f(b_k, shift)

				if `estimator'==2 {

					// shift identified from level ratio
					loc shiftres (_b[/g0]/_b[/b0]-1)

					forvalues k=0/`polynomial' {
						loc rhs (_b[/b`k']/(1+`shiftres'))
						loc teststr `teststr' (_b[/g`k'] = `rhs')
					}

				}
				else if `estimator'==3 {

					if "`log'"=="log" & `polynomial'>0 {
						loc shiftres ///
							(exp((_b[/g`=`polynomial'-1']-_b[/b`=`polynomial'-1']) ///
							/(`polynomial'*_b[/b`polynomial']))-1)
					}
					else {
						loc shiftres (_b[/g0]/_b[/b0]-1)
					}

					forvalues k=0/`polynomial' {

						if "`log'"=="" {

							// g_k = b_k (1+shift)^(k+1)

							loc rhs (_b[/b`k']*(1+`shiftres')^`=`k'+1')

						}
						else {

							// log-income case:
							// g_k = sum_{n=k}^K b_n comb(n,k) ln(1+shift)^(n-k)

							loc rhs _b[/b`k']

							if `polynomial'>`k' {
								forvalues n=`=`k'+1'/`polynomial' {
									loc rhs ///
									(`rhs' + ///
									_b[/b`n']*comb(`n',`k')*ln(1+`shiftres')^`=`n'-`k'')
								}
							}
						}

						loc teststr `teststr' (_b[/g`k'] = `rhs')
					}
				}
			}
			else if `estimator'==1 {

				if "`test'"!="notest" {
					forvalues k=0/`polynomial' {
						loc teststr `teststr' (_b[/b`k'] = _b[/g`k'])
					}
				}
			}
							
				
				//ESTIMATION AND INFERENCE
				tempname b V bu Vu tmpb
				
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
					
					//estimate unrestricted model
					`noisily' reg `unresmodel', nocons
					if `bootreps'==1 { //get variance for unrestricted model, perform test of model restrictions
						varcorrect `unresmodel', `smallsample'
						mat `Vu'=r(V)
						mat `bu'=e(b)
						mat colnames `bu'=`newnames'	
						mat colnames `Vu'=`newnames'	
						mat rownames `Vu'=`newnames'	
						ereturn post `bu' `Vu'
						mat `bu'=e(b)
						mat `Vu'=e(V)
						if "`test'"!="notest"&`estimator'>0 {
							`noisily' eret di
							`noisily' testnl `teststr'
							loc chi2=r(chi2)
							loc p_mod=r(p)
							loc df=r(df)
						}
					if "`saveunres'"!="" est save `saveunres', replace
					}
					else if "`transform'"!="notransform"|`estimator'>0 { //just change names and post
						mat `tmpb'=e(b)
						mat colnames `tmpb'=`newnames'
						ereturn post `tmpb'
					}
					if `s'>0 mat `Vu'=nullmat(`Vu') \ e(b)
					if `s'==0 mat `bu'=e(b)
					
					//ESTIMATE RESTRICTED MODEL 
					if inlist(`estimator',2,3) { //USING THE PROFILING ESTIMATOR
						profile23, ///
							y(`y') z(`z') bunch(`bunch') relbin(`relbin') ///
							cutoff_orig(`cutoff_orig') bw_orig(`bw_orig') //
							polynomial(`polynomial') estimator(`estimator') ///
							l(`L') h(`H') ///
							`log' `normalize
					}
					else if `estimator'==1 { //With OLS
						`noisily' reg `y' `rhsvars' `cons' `bunchvars', nocons
						mat `b'=e(b)
						mat colnames `b'=`namesest1'
						if `bootreps'==1 {
							varcorrect `y' `rhsvars' `cons' `bunchvars', `smallsample'
							mat `V'=r(V)
							mat rownames `V'=`namesest1'
							mat colnames `V'=`namesest1'
							ereturn post `b' `V'
							}
						else ereturn post `b'
					}
					
					////TRANSFORM ESTIMATES
					
					if "`transform'"!="notransform" {
						
						
						inlist(`estimator',1,0) { //TRANSFORM ESTIMATES
							local nlcomopt
							if "`nlcom'" != "" local nlcomopt nlcom(`nlcom')

							if inlist(`estimator', 2, 3) {
								if `s' == 0 {
									mat `b' = e(b)
									if `bootreps' == 1 mat `V' = e(V)
								}
								else if `bootreps' == 2 {
									mat `V' = nullmat(`V') \ e(b)
								}
							}
							else if "`transform'" != "notransform" {
								if `bootreps' == 1 {
									bunchcalc, estimator(`estimator') polynomial(`polynomial') ///
										cutoff(`cutoff_orig') bw(`bw_orig') h(`H') l(`L') ///
										b0(`b0') t0(`t0') t1(`t1') ///
										`constant' `positiveshift' `log' `nlcomopt' `normalize'
									mat `V' = r(V)
								}
								else {
									bunchcalc, estimator(`estimator') polynomial(`polynomial') ///
										cutoff(`cutoff_orig') bw(`bw_orig') h(`H') l(`L') ///
										b0(`b0') t0(`t0') t1(`t1') boot ///
										`constant' `positiveshift' `log' `nlcomopt' `normalize'
								}

								if r(exit) == 1 {
									noi di as error `"Could not find solution to polynomial equation for the response of the marginal buncher in one or more bootstrap repetitions. Try constant or notransform."'
									exit 301
								}

								if `s' == 0 {
									local nlcom `r(nlcom)'
									mat `b' = r(b)
								}
								else {
									mat `V' = nullmat(`V') \ r(b)
								}
							}
							else {
								if `s' == 0 {
									mat `b' = e(b)
									if `bootreps' == 1 mat `V' = e(V)
								}
								else if `bootreps' == 2 {
									mat `V' = nullmat(`V') \ e(b)
								}
							}

							if `s' > 0 noi _dots `s' 0
					}
					
					//bootstrap inference sunmmary & test
					if `bootreps'>1 {
						clear
						svmat `V'
						corr _all, cov
						mat `V'=r(C)
						if `estimator'>0 {
							clear 
							svmat `Vu'
							corr _all, cov
							mat `Vu'=r(C)
						}
					}
					
					//TEST RESTRICTIONS
					if `estimator'>0 {
						if "`test'"!="notest"&`bootreps'>0 {
							mat colnames `bu'=`newnames'
							mat colnames `Vu'=`newnames'
							mat rownames `Vu'=`newnames'
							ereturn post `bu' `Vu'		
							testnl `teststr'
							loc chi2=r(chi2)
							loc p_mod=r(p)
							loc df=r(df)
							}
						if "`saveunres'"!="" {
							mat colnames `bu'=subinstr("`names0'","delta","",.)
							if "`transform'"!="notransform" mat coleq `b'=`coleq0' `coleq1' bunching
							if `bootreps'>0 {
								mat colnames `Vu'=subinstr("`names0'","delta","",.)
								mat rownames `Vu'=subinstr("`names0'","delta","",.)
								if "`transform'"!="notransform" {
									mat coleq `V'=`coleq0' `coleq1' bunching
									mat roweq `V'=`coleq0' `coleq1' bunching	
								}
							ereturn post `bu' `Vu'
							}
							else ereturn post `bu'
							ereturn local cmd="regress"
							est sto `saveunres'
						}
					}

				
				//POST RESULTS
				mat colnames `b'=`names'
				if `bootreps'>=1 {
					mat colnames `V'=`names'
					mat rownames `V'=`names'
				}
				if "`transform'"!="notransform" {
					mat coleq `b'=`coleq0' `coleq1' bunching
					if `bootreps'>=1 {
						mat coleq `V'=`coleq0' `coleq1' bunching
						mat roweq `V'=`coleq0' `coleq1' bunching
					}
				}
				
				
				restore
				//return results
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
				ereturn scalar lower_limit=`zL'
				ereturn scalar upper_limit=`zH'
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
    syntax anything, [nl nosmallsample]

    quietly {
        gettoken y anything : anything

        tempvar res rss
        tempname g V
        local numbins = _N

        summarize `y', meanonly
        local N = r(sum)
        local k = e(df_m)

        if "`nl'" == "nl" {
            predictnl `res' = predict(), g(`g') iterate(1000)
            replace `res' = -`res'

            unab gvars : `g'*
            return local xvars "`gvars'"

            mata: st_matrix("`V'", ///
                varcorrect(st_data(., tokens(st_local("gvars"))), ///
                           st_data(., "`y'"), ///
                           st_data(., "`res'"), 0))
        }
        else {
            predict `res', residuals
            return local xvars "`anything'"

            mata: st_matrix("`V'", ///
                varcorrect(st_data(., tokens(st_local("anything"))), ///
                           st_data(., "`y'"), ///
                           st_data(., "`res'"), 0))
        }

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

program define profile23, eclass
    syntax, y(name) z(name) relbin(name) bunch(name) ///
        cutoff_orig(real) bw_orig(real) ///
        polynomial(integer) estimator(integer) l(integer) h(integer) ///
        [log nonormalize novar]

    mata: profile23_run( ///
        "`y'", "`z'", "`relbin'", "`bunch'", ///
        `cutoff_orig', `bw_orig', ///
        `polynomial', `estimator', ///
        ("`nonormalize'" == ""), ///
        ("`log'" == "log"), ///
        `l', `h', ///
        ("`var'" != "novar") ///
    )
	
	if "`var'"!="novar" ereturn post r_b_profile23 r_V_profile23
	else  ereturn post r_b_profile23
end

mata:

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

        e_i = e_stack
        e_i[s] = e_i[s] + N

        meat = meat + fw_orig[i] * (G' * e_i * e_i' * G)
    }

    bread = invsym(quadcross(G, G) :* N)

    return(bread * meat * bread)
}

void profile23_run(
    string scalar yvar,
    string scalar zvar,
    string scalar relbinvar,
    string scalar bunchvar,
    real scalar cutoff_orig,
    real scalar bw_orig,
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar L,
    real scalar H,
    real scalar dovar
)
{
    real colvector y, z, relbin, bunch, side
    real colvector yout, zout, sideout, ystack, resid
    real colvector stack_id, fw_orig
    real scalar Hstar, lndelta_hat, delta_hat
    real scalar i, idx, nleft, nright
    real rowvector beta_hat, b
    real matrix X, dX, G, Vout
    transmorphic S

    y      = st_data(., yvar)
    z      = st_data(., zvar)
    relbin = st_data(., relbinvar)
    bunch  = st_data(., bunchvar)

    side = J(rows(y), 1, .)
    side[relbin :< -L + 1] = -1
    side[relbin :>  H]     =  1

    Hstar = sum(select(y, bunch :> 0))

    S = optimize_init()
    optimize_init_evaluator(S, &profQ23())
    optimize_init_evaluatortype(S, "d0")
    optimize_init_params(S, ln(0.05))

    optimize_init_argument(S, 1,  y)
    optimize_init_argument(S, 2,  z)
    optimize_init_argument(S, 3,  side)
    optimize_init_argument(S, 4,  bunch)
    optimize_init_argument(S, 5,  Hstar)
    optimize_init_argument(S, 6,  cutoff_orig)
    optimize_init_argument(S, 7,  bw_orig)
    optimize_init_argument(S, 8,  K)
    optimize_init_argument(S, 9,  estimator)
    optimize_init_argument(S, 10, normalized)
    optimize_init_argument(S, 11, islog)
    optimize_init_argument(S, 12, L)
    optimize_init_argument(S, 13, H)

    lndelta_hat = optimize(S)
    delta_hat   = exp(lndelta_hat)

    beta_hat = profBeta23(
        delta_hat,
        y, z, side, bunch, Hstar,
        cutoff_orig, bw_orig,
        K, estimator, normalized, islog, L, H
    )

    b = beta_hat, delta_hat, Hstar

    if (dovar == 1) {
        yout    = select(y,    bunch :== 0)
        zout    = select(z,    bunch :== 0)
        sideout = select(side, bunch :== 0)

        X = makeX23(delta_hat, zout, sideout, cutoff_orig, bw_orig, K, ///
                    estimator, normalized, islog, L, H)

        ystack = select(yout, sideout :== -1) \ ///
                 select(yout, sideout :==  1) \ ///
                 Hstar

        resid = ystack - X * beta_hat'

        dX = dXddelta23(delta_hat, zout, sideout, cutoff_orig, bw_orig, K, ///
                        estimator, normalized, islog, L, H)

        // Jacobian wrt beta, delta, Hstar
        G = X, dX * beta_hat', J(rows(X), 1, 0)
        G[rows(G), cols(G)] = 1

        // Map original bins to stacked rows:
        // outside-left bins first, outside-right bins second, excluded bins last.
        nleft  = sum((bunch :== 0) :& (side :== -1))
        nright = sum((bunch :== 0) :& (side :==  1))

        stack_id = J(rows(y), 1, .)
        idx = 0

        for (i = 1; i <= rows(y); i++) {
            if (bunch[i] == 0 & side[i] == -1) {
                idx = idx + 1
                stack_id[i] = idx
            }
        }

        for (i = 1; i <= rows(y); i++) {
            if (bunch[i] == 0 & side[i] == 1) {
                idx = idx + 1
                stack_id[i] = idx
            }
        }

        for (i = 1; i <= rows(y); i++) {
            if (bunch[i] > 0) {
                stack_id[i] = nleft + nright + 1
            }
        }

        fw_orig = y

        Vout = varcorrect_collapsed(G, fw_orig, stack_id, resid, 0)
    }
    else {
        Vout = J(cols(b), cols(b), .)
    }

    st_matrix("r_b_profile23", b)
    st_matrix("r_V_profile23", Vout)
}

real scalar profQ23(
    real scalar lndelta,
    real colvector y,
    real colvector z,
    real colvector side,
    real colvector bunch,
    real scalar Hstar,
    real scalar cutoff_orig,
    real scalar bw_orig,
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar L,
    real scalar H
)
{
    real scalar delta
    real colvector yout, zout, sideout, ystack, resid
    real matrix X, W, beta

    delta = exp(lndelta)

    yout    = select(y,     bunch :== 0)
    zout    = select(z,     bunch :== 0)
    sideout = select(side,  bunch :== 0)

    X = makeX23(delta, zout, sideout, cutoff_orig, bw_orig, K, ///
                estimator, normalized, islog, L, H)

    ystack = select(yout, sideout :== -1) \ ///
             select(yout, sideout :==  1) \ ///
             Hstar

    beta =  qrsolve(X, ystack)
    resid = ystack - X * beta

    return(quadcross(resid,resid))
}


real rowvector profBeta23(
    real scalar delta,
    real colvector y,
    real colvector z,
    real colvector side,
    real colvector bunch,
    real scalar Hstar,
    real scalar cutoff_orig,
    real scalar bw_orig,
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar L,
    real scalar H
)
{
    real colvector yout, zout, sideout, ystack
    real matrix X, beta

    yout    = select(y,     bunch :== 0)
    zout    = select(z,     bunch :== 0)
    sideout = select(side,  bunch :== 0)

    X = makeX23(delta, zout, sideout, cutoff_orig, bw_orig, K, ///
                estimator, normalized, islog, L, H)

    ystack = select(yout, sideout :== -1) \ ///
             select(yout, sideout :==  1) \ ///
             Hstar

    beta =  qrsolve(X, ystack)

    return(beta')
}

real matrix dXddelta23(
    real scalar delta,
    real colvector z,
    real colvector side,
    real scalar cutoff_orig,
    real scalar bw_orig,
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar L,
    real scalar H
)
{
    real scalar eps, dm, dp
    real matrix Xp, Xm

    eps = max((1e-6, abs(delta)*1e-5))
    dp  = delta + eps
    dm  = max((delta - eps, 1e-10))

    Xp = makeX23(dp, z, side, cutoff_orig, bw_orig, K, ///
                 estimator, normalized, islog, L, H)

    Xm = makeX23(dm, z, side, cutoff_orig, bw_orig, K, ///
                 estimator, normalized, islog, L, H)

    return((Xp - Xm)/(dp - dm))
}

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


real colvector transform_right_z(
    real colvector zR,
    real scalar delta,
    real scalar cutoff_orig,
    real scalar bw_orig,
    real scalar normalized,
    real scalar islog
)
{
    if (islog == 1) {
        return(zR :- ln(1 + delta))
    }

    if (normalized == 1) {
        return(((cutoff_orig :+ bw_orig*zR) :/ (1+delta) :- cutoff_orig) :/ bw_orig)
    }

    return(zR :/ (1 + delta))
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


real matrix makeX23(
    real scalar delta,
    real colvector z,
    real colvector side,
    real scalar cutoff_orig,
    real scalar bw_orig,
    real scalar K,
    real scalar estimator,
    real scalar normalized,
    real scalar islog,
    real scalar L,
    real scalar H
)
{
    real scalar ex_lo, ex_hi, resp, cutoff_est, zmax_est
    real colvector zL, zR, zR0
    real matrix X0, X1, XB, X

    zL = select(z, side :== -1)
    zR = select(z, side :==  1)

    // Left side: h0(z)
    X0 = pbasis(zL, K)

    // Right side
    if (estimator == 2) {

        // Chetty-style coefficient restriction:
        // g_k = b_k/(1+delta)

        X1 = pbasis(zR, K) / (1 + delta)

    }
    else {

        // Theory-consistent restriction:
        // h1(z)=h0(z/(1+delta))/(1+delta)

        zR0 = transform_right_z(
            zR,
            delta,
            cutoff_orig,
            bw_orig,
            normalized,
            islog
        )

        X1 = pbasis(zR0, K) / (1 + delta)
    }

    // Estimation-scale limits
    if (normalized == 1) {

        cutoff_est = 0

        ex_lo = -L + 1
        ex_hi =  H + 1

        // matches old:
        // zmax = max(z)+bw/2
        // with normalized bw=1
        zmax_est = max(z) + 0.5
    }
    else {

        cutoff_est = cutoff_orig

        ex_lo = cutoff_orig + (-L + 1)*bw_orig
        ex_hi = cutoff_orig + ( H + 1)*bw_orig

        zmax_est = max(z) + bw_orig/2
    }

    // Response length
    resp = response_length(
        delta,
        cutoff_orig,
        bw_orig,
        normalized,
        islog
    )

    if (estimator == 2) {

        /*
        Profile analogue of old NLS estimator 2:

            B*bw = delta * ∫_{cutoff}^{zmax} h0(z) dz

        Therefore

            H* =
                ∫excluded h0(z) dz
                +
                delta * ∫_{cutoff}^{zmax} h0(z) dz
        */

        XB =
            intbasis(ex_lo, ex_hi, K)
            +
            delta * intbasis(cutoff_est, zmax_est, K)

    }
    else {

        /*
        Estimator 3:

            H* =
                ∫excluded h0(z) dz
                +
                ∫_{cutoff}^{cutoff+response} h0(z) dz
        */

        XB =
            intbasis(ex_lo, ex_hi, K)
            +
            intbasis(cutoff_est,
                     cutoff_est + resp,
                     K)
    }

    X = X0 \ X1 \ XB

    return(X)
}

end
