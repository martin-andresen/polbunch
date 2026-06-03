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
								noi di in red "Could not estiamte separate polynomials on either side of the cutoff."
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
					foreach k of numlist 1/`polynomial' 0 {
						loc newnames `newnames' /`l'`k'
						if `estimator'==1&"`l'"=="b" loc namesest1 `namesest1' /b`k'
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

				//NAMES for final model
				if "`transform'"!="notransform" {
					if `estimator'!=2 loc names `names' `names' number_bunchers excess_mass shift marginal_response
					else loc names `names' `names' delta number_bunchers excess_mass shift marginal_response
					if "`t0'"!=""&"`t1'"!="" loc names `names' elasticity
				}
				else {
					tokenize `names'
					loc names
					forvalues i=1/`=`L'' {
						loc bnames `bnames' number_bunchers:`=`L'-`i'+1'.below
					}
					if `H'>0 {
							forvalues i=1/`=`H'' {
								loc bnames `bnames' number_bunchers:`i'.above
							}
						}
					if `estimator'>1 {
						if `estimator'==2 {
							if "`positiveshift'"!="nopositiveshift" loc names `names' delta:lndelta
							else loc names `names' delta:delta
						}
						else {
							if "`positiveshift'"!="nopositiveshift" loc names `names' shift:lnshift
							else loc names `names' shift:shift
						}
						forvalues k=1/`polynomial' {
							loc names `names' h0:``k''
						}
						loc names `bnames' `names'
					} 
				else {
					forvalues h=0/1 {
						forvalues k=1/`polynomial' {
							loc names `names' h`h':``k''
						}
						loc names `names' h`h':_cons
						if `estimator'==1 continue, break
					}
				loc names `names' `bnames'
				}
				
				}
				
				
				
				if `estimator'>1 { //specify NL model
					
				//Specify parameter restrictions
				
				*********shift*********
				if "`positiveshift'"!="nopositiveshift" loc shift exp({lnshift})
				else loc shift {shift}
				
				*********b0************
				loc B=0
				forvalues b=1/`=`H'+`L'' {
					loc B `B'+{bunch`b'}
					}
				loc cns (`B')*`bw'
				if `estimator'!=2 {
					forvalues k=1/`polynomial' {
						if "`log'"=="" loc cns `cns' - ((`cutoff'^`=`k'+1')*((1+`shift')^(`=`k'+1')-1)*({b`k'}/`=`k'+1'))
						else loc cns `cns' - (((`cutoff'+ln(1+`shift'))^(`=`k'+1')-`cutoff'^`=`k'+1')*({b`k'}/`=`k'+1'))
						}
						
					if "`log'"=="" loc b0 (`cns')/(`cutoff'*`shift')
					else loc b0 (`cns')/ln(1+`shift')
				}
				else {
					loc cns ((`cns')/`shift')
					su `z'
					loc zmax=r(max)+`bw'/2				
					forvalues k=1/`polynomial' {
						loc cns `cns' - ({b`k'}/`=`k'+1')*(`zmax'^`=`k'+1'-`cutoff'^`=`k'+1')
					}
				loc b0 ((`cns')/(`zmax'-`cutoff'))
				}
				
				********b1-bK*********
				forvalues k=1/`polynomial' {
					loc b`k' {b`k'}
				}
				
				********g0-gK*********
				forvalues k=0/`polynomial' {
					if `estimator'==1 loc g`k' `b`k''
					else if `estimator'==2 loc g`k' (`b`k'')/(1+`shift')
					else if "`log'"=="" loc g`k' (`b`k'')*(1+`shift')^`=`k'+1'
					else {
						loc g`k' `b`k''
						if `polynomial'>`k' forvalues n=`=`k'+1'/`polynomial' {
							loc g`k' `g`k'' +{b`n'}*comb(`n',`k')*ln(1+`shift')^`=`n'-`k''
						}
					}
				}
				
					
				**** NL model string ****
				loc modstr (`b0')*0.`dum'
				forvalues k=1/`polynomial' {
					loc modstr `modstr' +(`b`k'')*0.`dum'*`z'^`k'
				}
				loc modstr `modstr' +(`g0')*1.`dum'
				forvalues k=1/`polynomial' {
					loc modstr `modstr' +(`g`k'')*1.`dum'*`z'^`k'
				}
				
				**** bunch dummies ****
				forvalues bval=1/`=`H'+`L'' {
					loc modstr `modstr' + {bunch`bval'}*`bval'.`bunch'
					}

				
				//REPLACE { } with _b[/ ] for transformations and testing
				foreach l in b g {
					forvalues k=0/`polynomial' {
						loc `l'`k' = subinstr("``l'`k''","{","_b[/",.)
						loc `l'`k' = subinstr("``l'`k''","}","]",.)
					}
				}

				//TESTSTRING
				if "`test'"!="notest" {
						if `estimator'==3&"`log'"=="log" loc shiftres (exp((_b[/g`=`polynomial'-1']-_b[/b`=`polynomial'-1'])/(`polynomial'*_b[/b`polynomial']))-1)
						else loc shiftres (_b[/g0]/_b[/b0]-1)
						forvalues k=0/`polynomial' {
							if "`positiveshift'"=="nopositiveshift" loc teststr `teststr' (_b[/g`k']=`=subinstr("`g`k''","_b[/shift]","`shiftres'",.)')
							else loc teststr `teststr' (_b[/g`k']=`=subinstr("`g`k''","exp(_b[/lnshift])","`shiftres'",.)')
						}
					}
				}
				else {
					loc b0 _b[/b0]
					if `estimator'==1&"`test'"!="notest" {
						forvalues k=0/`polynomial' {
							loc teststr `teststr' (_b[/b`k']=_b[g`k'])
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
					if `estimator'>1 { //with NLS
						// find good starting values
						// Use excess mass relative to predicted left density at cutoff.

						local Binit = 0
						forvalues bval = 1/`=`H'+`L'' {
							local thisb = _b[/bunch`bval']
							if !missing(`thisb') {
								local Binit = `Binit' + max(0, `thisb')
							}
						}

						// fallback-safe predicted left-side bin count at cutoff
						capture local hinit = _b[/b0]
						if _rc {
							quietly summarize `y' if `z' < `cutoff' & `bunch' == 0, meanonly
							local hinit = r(mean)
						}

						forvalues k = 1/`polynomial' {
							local hinit = `hinit' + _b[/b`k']*(`cutoff'^`k')
						}

						// fallback if polynomial prediction is bad
						if missing(`hinit') | `hinit' <= 0 {
							quietly summarize `y' if `z' < `cutoff' & `bunch' == 0, meanonly
							local hinit = r(mean)
						}

						// approximate shift: excess bins / local counterfactual bin height
						local shiftinit = (`Binit' * `bw') / (`hinit' * `cutoff')

						if missing(`shiftinit') | `shiftinit' <= 0 {
							local shiftinit = 0.05
						}
						if `shiftinit' > 0.5 {
							local shiftinit = 0.5
						}

						// initialize shift
						if "`positiveshift'" != "nopositiveshift" {
							local init lnshift `=ln(`shiftinit')'
						}
						else {
							local init shift `shiftinit'
						}

						// initialize b0; needed especially for polynomial(0)
						capture local b0init = _b[/b0]
						if _rc {
							quietly summarize `y' if `z' < `cutoff' & `bunch' == 0, meanonly
							local b0init = r(mean)
						}
						local init `init' b0 `b0init'

						// initialize higher-order polynomial terms; skipped automatically if polynomial==0
						forvalues k = 1/`polynomial' {
							local init `init' b`k' `=_b[/b`k']'
						}

						// initialize bunching dummy parameters
						forvalues bval = 1/`=`H'+`L'' {
							local init `init' bunch`bval' `=_b[/bunch`bval']'
						}
						
						
						//estimate restricted model W/ NLS
						`noisily' nl (`y'=`modstr'), init(`init')
						if `bootreps'==1 {
							varcorrect `y', nl `smallsample'
							
							mat `V'=r(V)
							ereturn repost V=`V'
							}
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
						if "`transform'"!="notransform" { //TRANSFORM ESTIMATES
							if `bootreps'==1 {
								bunchcalc, estimator(`estimator') polynomial(`polynomial') cutoff(`cutoff') bw(`bw') h(`H') l(`L') b0(`b0') t0(`t0') t1(`t1') `constant' `positiveshift' `log' nlcom(`nlcom')
								mat `V'=r(V)
							}
							else bunchcalc, estimator(`estimator') polynomial(`polynomial') cutoff(`cutoff') bw(`bw') h(`H') l(`L') b0(`b0') t0(`t0') t1(`t1') boot `constant' `positiveshift' `log' nlcom(`nlcom')
							if r(exit)==1 {
								noi di in red "Could not find solution to polynomial equation for the response of the marginal buncher in one or more bootstrap repetitions. You could consider trying the constant approximation using the option "constant", or the option "notransform" to report raw estimates and then manually convert those to objects of interest post-estimation (the latter being less prone to bias)."
								exit 301
							}
							if `s'==0 {
								loc nlcom `=r(nlcom)'
								mat `b'=r(b)
							}
							else {
								mat `V'=nullmat(`V') \ r(b)
							}
						}
						else if `s'==0 {
							mat `b'=e(b)
							if `bootreps'==1 mat `V'=e(V)
							}
						else if `s'>0&`bootreps'==2 mat `V'=nullmat(`V') \ e(b)
						if `s'>0 noi _dots `s' 0
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


real scalar eresp(real scalar B, real scalar tau, real matrix cf, real scalar bw)
{
    real matrix integral, roots, realroots, out

    if (cols(cf) == 1) {
        if (cf[1] <= 0) return(.)
        return((B * bw) / cf[1])
    }

    integral = polyinteg(cf, 1)
    integral[1] = -polyeval(integral, tau) - B*bw

    roots = polyroots(integral)
    realroots = Re(select(roots, Im(roots) :== 0))
    out = sort(select(realroots, realroots :> tau)', 1)'

    if (cols(out) == 0) return(.)
    else return(out[1] - tau)
}

end
end

