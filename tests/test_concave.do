* test_concave.do -- concave-kink (t1<t0, bunchers from BELOW the cutoff) support
*
* Run from the repository root:   do tests/test_concave.do
*
*   A. convex output unchanged (regression values)
*   B. log mode: the concave fit is the exact mirror image of the convex fit
*      (coefficients, standard errors, elasticity, polbunchbias outputs)
*   C. level mode: estimator 3 recovers the true elasticity; the inline bias
*      correction recovers it for estimators 1 and 5 (Monte Carlo)
*   D. downstream commands run on a concave fit (plot, contrast, permute)
clear all
adopath ++ "`c(pwd)'"

scalar nfail = 0
program define _ck
    gettoken msg rest : 0
    local ok = `rest'
    if `ok' di as txt "  ok    `msg'"
    else {
        di as error "  FAIL  `msg'"
        scalar nfail = nfail + 1
    }
end

*---------------------------------------------------------------- A. convex
di as txt _n "A. convex orientation unchanged"
set seed 4242
polbunchgendata z, obs(40000) cutoff(1) el(0.4) t0(0.2) t1(0.6) distribution(0.2+rbeta(2,3)*3)
qui polbunch z, cutoff(1) bw(0.03) poly(3) estimator(3) t0(0.2) t1(0.6) window(3 3) vce(none)
_ck "convex est3 elasticity = .399163928757" (abs(_b[bunching:elasticity] - .399163928757) < 1e-9)
qui polbunch z, cutoff(1) bw(0.03) poly(3) estimator(1) t0(0.2) t1(0.6) window(3 3) vce(none)
_ck "convex est1 bias_elasticity = .0505381188808" (abs(e(bias_elasticity) - .0505381188808) < 1e-9)
_ck "e(concave) = 0" (e(concave) == 0)

*---------------------------------------------------------------- B. log mirror
di as txt _n "B. log mode: concave fit == mirror image of the convex fit"
foreach est in 1 3 5 {
    clear
    set seed 31
    polbunchgendata z, obs(100000) cutoff(0) el(0.4) t0(0.2) t1(0.6) distribution(rbeta(2,2)*4-2) log
    qui polbunch z, cutoff(0) bw(0.03) poly(3) estimator(`est') t0(0.2) t1(0.6) log window(3 12) nobias
    local e1 = _b[bunching:elasticity]
    local s1 = _se[bunching:elasticity]
    local m1 = _b[bunching:excess_mass]
    local h1 = _b[h0:z]
    local mr1 = _b[bunching:marginal_response]
    if inlist(`est', 1, 5) {
        qui polbunchbias
        local bb1 = r(bias_elasticity)
        local bs1 = r(bias_slope)
    }
    replace z = -z
    qui polbunch z, cutoff(0) bw(0.03) poly(3) estimator(`est') t0(0.6) t1(0.2) log window(12 3) nobias
    _ck "est`est' concave flag"                (e(concave) == 1)
    _ck "est`est' elasticity identical"        (abs(_b[bunching:elasticity] - `e1') < 1e-9)
    _ck "est`est' se(elasticity) identical"    (abs(_se[bunching:elasticity] - `s1') < 1e-9)
    _ck "est`est' excess_mass identical"       (abs(_b[bunching:excess_mass] - `m1') < 1e-9)
    _ck "est`est' odd coefficient flips sign"  (abs(_b[h0:z] + `h1') < 1e-9)
    _ck "est`est' marginal_response flips"     (abs(_b[bunching:marginal_response] + `mr1') < 1e-9)
    _ck "est`est' e(t0)/e(t1) are the user's"  (e(t0) == 0.6 & e(t1) == 0.2)
    _ck "est`est' excluded region mirrored"    (abs(e(lower_limit) + 0.36) < 1e-9 & abs(e(upper_limit) - 0.09) < 1e-9)
    if inlist(`est', 1, 5) {
        qui polbunchbias
        _ck "est`est' polbunchbias bias_elasticity identical" (abs(r(bias_elasticity) - `bb1') < 1e-9)
        _ck "est`est' polbunchbias bias_slope flips"          (abs(r(bias_slope) + `bs1') < 1e-9)
        _ck "est`est' polbunchbias t0/t1 restored"            (r(t0) == 0.6 & r(t1) == 0.2)
    }
}

*---------------------------------------------------------------- C. MC truth
di as txt _n "C. Monte Carlo against a concave DGP (true elasticity 0.4)"
foreach lg in "" "log" {
    local n1 = 0
    local s3 = 0
    local sc1 = 0
    local sc5 = 0
    forvalues r = 1/12 {
        clear
        set seed `=2000+`r''
        if "`lg'" == "" polbunchgendata z, obs(80000) cutoff(1) el(0.4) t0(0.6) t1(0.2) distribution(0.2+rbeta(2,3)*3)
        else polbunchgendata z, obs(80000) cutoff(0) el(0.4) t0(0.6) t1(0.2) distribution(rbeta(2,2)*4-2) log
        local cut = cond("`lg'" == "log", 0, 1)
        local ++n1
        qui polbunch z, cutoff(`cut') bw(0.03) poly(3) estimator(3) t0(0.6) t1(0.2) window(10 3) vce(none) nobias `lg'
        local s3 = `s3' + _b[bunching:elasticity]
        qui polbunch z, cutoff(`cut') bw(0.03) poly(3) estimator(1) t0(0.6) t1(0.2) window(10 3) vce(none) `lg'
        local sc1 = `sc1' + _b[bunching:elasticity] - e(bias_elasticity)
        qui polbunch z, cutoff(`cut') bw(0.03) poly(3) estimator(5) t0(0.6) t1(0.2) window(10 3) vce(none) `lg'
        local sc5 = `sc5' + _b[bunching:elasticity] - e(bias_elasticity)
    }
    local lab = cond("`lg'" == "log", "log  ", "level")
    _ck "`lab' est3 mean elasticity ~ 0.4 (`=string(`s3'/`n1',"%6.4f")')"            (abs(`s3'/`n1' - 0.4) < 0.01)
    _ck "`lab' est1 bias-corrected ~ 0.4 (`=string(`sc1'/`n1',"%6.4f")')"            (abs(`sc1'/`n1' - 0.4) < 0.015)
    _ck "`lab' est5 bias-corrected ~ 0.4 (`=string(`sc5'/`n1',"%6.4f")')"            (abs(`sc5'/`n1' - 0.4) < 0.015)
}

*---------------------------------------------------------------- D. downstream
di as txt _n "D. downstream commands on a concave fit"
clear
set seed 777
polbunchgendata z, obs(100000) cutoff(1) el(0.4) t0(0.6) t1(0.2) distribution(0.2+rbeta(2,3)*3)
qui polbunch z, cutoff(1) bw(0.03) poly(3) estimator(3) t0(0.6) t1(0.2) window(10 3) nobias
estimates store m3
qui polbunch z, cutoff(1) bw(0.03) poly(5) estimator(1) t0(0.6) t1(0.2) window(10 3) nobias
estimates store m1
capture polbunchplot m3
_ck "polbunchplot (single)" (_rc == 0)
capture polbunchplot m3 m1
_ck "polbunchplot (two models)" (_rc == 0)
estimates restore m3
capture polbunchplot, rootogram
_ck "polbunchplot, rootogram" (_rc == 0)
capture polbunch_contrast m3 m1, reps(20) nodots
_ck "polbunch_contrast" (_rc == 0)
* permutation test needs room around the fitting window
clear
set seed 777
polbunchgendata z, obs(300000) cutoff(3) el(0.3) t0(0.6) t1(0.2) distribution(3+rnormal(0,1))
qui polbunch z if inrange(z,2,4), cutoff(3) bw(0.03) poly(3) estimator(3) t0(0.6) t1(0.2) window(8 1) nobias
capture polbunch_permute, maxcutoffs(15) nodots
_ck "polbunch_permute" (_rc == 0)

* savebins() + pre-binned vce(cluster stub): round trip in the user's orientation
clear
set seed 5
polbunchgendata z, obs(60000) cutoff(1) el(0.4) t0(0.6) t1(0.2) distribution(0.2+rbeta(2,3)*3)
gen id = ceil(runiform()*400)
qui polbunch z, cutoff(1) bw(0.03) poly(3) estimator(1) t0(0.6) t1(0.2) window(10 3) vce(cluster id) nobias savebins("`c(tmpdir)'/pb_sb_concave", replace)
scalar rt_el = _b[bunching:elasticity]
scalar rt_se = _se[bunching:elasticity]
use "`c(tmpdir)'/pb_sb_concave", clear
qui polbunch freq midpoint, cutoff(1) poly(3) estimator(1) t0(0.6) t1(0.2) window(10 3) vce(cluster cv) nobias
_ck "savebins round trip: elasticity" (abs(_b[bunching:elasticity]-rt_el) < 1e-6)
_ck "savebins round trip: se" (abs(_se[bunching:elasticity]-rt_se) < 1e-6)

di as txt _n "failures: " nfail
if nfail > 0 exit 9
