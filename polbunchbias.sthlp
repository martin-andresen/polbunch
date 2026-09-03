{smcl}
{* *! version 2.0.0 26sep2026}{...}
{vieweralsosee "[R] return" "help return"}{...}
{viewerjumpto "Syntax" "polbunchbias##syntax"}{...}
{viewerjumpto "Description" "polbunchbias##description"}{...}
{viewerjumpto "Options" "polbunchbias##options"}{...}
{viewerjumpto "Examples" "polbunchbias##examples"}{...}
{viewerjumpto "Stored results" "polbunchbias##results"}{...}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{hi:polbunchbias} {hline 2}}Analytical bias of polynomial bunching estimators{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Post-estimation mode (immediately after {cmd:polbunch}):

{p 8 15 2}
{cmd:polbunchbias}
[{it:options}]

{pstd}
Standalone mode (all required options supplied explicitly):

{p 8 15 2}
{cmd:polbunchbias}
{cmd:,}
{cmd:estimator(}{it:#}{cmd:)}
{cmd:zstar(}{it:#}{cmd:)}
{cmd:t0(}{it:#}{cmd:)}
{cmd:t1(}{it:#}{cmd:)}
{cmd:elasticity(}{it:#}{cmd:)}
{cmd:zlo(}{it:#}{cmd:)}
{cmd:zhi(}{it:#}{cmd:)}
{cmd:zl(}{it:#}{cmd:)}
{cmd:zh(}{it:#}{cmd:)}
{c -(}{cmd:relslope(}{it:#}{cmd:)} {cmd:|} {cmd:h0poly(}{it:numlist}{cmd:)}{c )-}
[{it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required (standalone mode only)}
{synopt :{cmd:estimator(}{it:#}{cmd:)}}estimator whose bias is evaluated: 0, 1, 2, 3 or 4{p_end}
{synopt :{cmd:zstar(}{it:#}{cmd:)}}bunching point or kink point, z*{p_end}
{synopt :{cmd:t0(}{it:#}{cmd:)}}lower tax rate, t0{p_end}
{synopt :{cmd:t1(}{it:#}{cmd:)}}higher tax rate, t1{p_end}
{synopt :{cmd:elasticity(}{it:#}{cmd:)}}true elasticity against which bias is evaluated{p_end}
{synopt :{cmd:zlo(}{it:#}{cmd:)}}lower bound of the estimation (fitting) window{p_end}
{synopt :{cmd:zhi(}{it:#}{cmd:)}}upper bound of the estimation (fitting) window{p_end}
{synopt :{cmd:zl(}{it:#}{cmd:)}}lower edge of the excluded / bunching region{p_end}
{synopt :{cmd:zh(}{it:#}{cmd:)}}upper edge of the excluded / bunching region{p_end}
{synopt :{cmd:relslope(}{it:#}{cmd:)}}relative slope of a degree-1 counterfactual at z*, h0'(z*)/h0(z*); mutually exclusive with {cmd:h0poly()}{p_end}
{synopt :{cmd:h0poly(}{it:numlist}{cmd:)}}counterfactual density as a polynomial in the centred running variable s, {it:lowest order first} (b0 = height at z*, b1 = slope, ...); any degree{p_end}

{pstd}
In {it:post-estimation mode} none of the above are given; all are read from
{cmd:e()}.  {cmd:relslope()} or {cmd:h0poly()} may still be supplied to
override the counterfactual read from {cmd:e(b)}.

{synoptset 24 tabbed}{...}
{syntab:Optional (both modes)}
{synopt :{cmd:log}}calculate in log z rather than level z (post-estimation mode reads this from {cmd:e(log)}){p_end}
{synopt :{cmd:constant}}invert bunching mass to the response with the constant-density shortcut (default for estimators 1 and 2){p_end}
{synopt :{cmd:exact}}invert bunching mass to the response by solving the density integral (default for estimators 0 and 3){p_end}
{synopt :{cmd:poolmass}}back out bunching mass by pooling the fitted counterfactual over the whole excluded window (Chetty 2011 eq. 16; default for estimators 1 and 2){p_end}
{synopt :{cmd:splitmass}}back out bunching mass by splitting at the kink (default for estimators 0 and 3){p_end}
{synopt :{cmd:bw(}{it:#}{cmd:)}}bin width used to scale bunching mass; default {cmd:e(bw_orig)} post-estimation, {cmd:bw(1)} standalone{p_end}
{synopt :{cmd:iterate}}solve the self-consistent fixed point theta_hat = theta_true + bias(theta_true){p_end}
{synopt :{cmd:tolerance(}{it:#}{cmd:)}}convergence tolerance for {cmd:iterate}; default {cmd:tolerance(1e-10)}{p_end}
{synopt :{cmd:maxiter(}{it:#}{cmd:)}}maximum iterations for {cmd:iterate}; default {cmd:maxiter(100)}{p_end}
{synopt :{cmd:underrelax(}{it:#}{cmd:)}}under-relaxation factor for {cmd:iterate}; default {cmd:underrelax(0.5)}{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:polbunchbias} computes the analytical bias of the polynomial bunching
estimators implemented by {helpb polbunch} (estimator codes 0, 1, 2, 3, 4).
The counterfactual density h0 is treated as a polynomial of arbitrary degree
K in the centred running variable {it:s} ({it:s} = z - z* in levels,
{it:s} = ln z - ln z* in logs).  {ul:When the true counterfactual is a
degree-K polynomial the reported biases are exact}, not a first-order
approximation.  Supplying a degree-1 polynomial (a constant and a slope)
reproduces the earlier local-linear behaviour.

{pstd}
The command reports the bias in each fitted h0 coefficient (b0, b1, ...,
bK) and in the downstream estimands: the relative slope, the estimated
bunching mass, the marginal response, the proportional shift, and the
elasticity.

{title:Modes}

{pstd}
{it:Post-estimation mode.}  Called immediately after {cmd:polbunch} with none
of the nine primary options, {cmd:polbunchbias} reads the kink point, bin
width, tax rates, elasticity, estimation window, excluded-region edges and
log/level flag from {cmd:e()}.  For the polynomial estimators (0-3) it also
reads the full fitted h0 polynomial from {cmd:e(b)} (degree
{cmd:e(polynomial)}) and recentres it on the cutoff, so the reported biases
are exact for the fitted order.  For the Saez estimator (4) the slope is
imputed from the two-point counterfactual in {cmd:_b[h0:_cons]} /
{cmd:_b[h1:_cons]}.  {cmd:relslope()} or {cmd:h0poly()} override the
counterfactual read from {cmd:e(b)}.

{pstd}
{it:Standalone mode.}  All nine primary options plus exactly one of
{cmd:relslope()} or {cmd:h0poly()} are supplied and only those values are
used.  A partial specification is an error.

{pstd}
{it:Log mode.}  With {cmd:log} (or {cmd:e(log)==1}), and exactly as in
{helpb polbunch}, the running variable is assumed to be {ul:already} in logs:
{cmd:zstar()}, {cmd:zlo/zhi()}, {cmd:zl/zh()} and the {cmd:h0poly()} /
{cmd:relslope()} coefficients are all in ln-earnings units.  The
counterfactual is then a polynomial in {it:s} = ln z - ln z*, the behavioural
response is the pure log translation {it:h1(s) = h0(s + rho)}, and there is
no proportional-shift bias.  Exact for any degree K, just like levels.

{title:The counterfactual and the response}

{pstd}
{cmd:h0poly(}{it:b0 b1 b2 ...}{cmd:)} gives h0 as a polynomial in the centred
running variable {it:s}, lowest order first: h0(s) = b0 + b1 s + b2 s^2 + ...
b0 is the height at the cutoff (any positive value; the response and
elasticity biases are scale-free).  {cmd:relslope(}{it:#}{cmd:)} is shorthand
for a degree-1 counterfactual with b0 = 1 and relative slope
{it:# = h0'/h0} at the cutoff, w.r.t. the running variable (per unit z in
levels, per unit ln z in logs).

{pstd}
The tax change implies {it:tau} = (1 - {it:t0})/(1 - {it:t1}); the true
response follows from {cmd:elasticity()}: in levels the proportional shift is
{it:Delta = tau^elasticity - 1}, in logs the log response is
{it:rho = ln(tau^elasticity)}.

{pstd}
For a polynomial estimator the degree K also fixes the polynomial order the
estimator is modelled as fitting -- these always coincide (post-estimation
mode reads {cmd:e(polynomial)}).  If the fitting window is too short to
identify a degree-K polynomial the fitting-stage biases are returned as
missing with a note.


{marker options}{...}
{title:Options}

{phang}
{cmd:estimator(}{it:#}{cmd:)} selects the estimator, matching the
{cmd:estimator()} codes of {helpb polbunch}:

{pmore}
{cmd:0} -- fits h0 from the left (uncontaminated) data only and h1 freely
from the right; {cmd:3} -- imposes the correct shift-and-stretch
restriction.  Both have no fitting-stage bias; they carry a bias only under
{cmd:constant} and/or {cmd:poolmass}.

{pmore}
{cmd:1} -- fits one counterfactual polynomial to the left and right included
regions and reads bunching mass as observed excluded-region mass minus the
fitted counterfactual.

{pmore}
{cmd:2} -- Chetty-style estimator restricting h1 = h0/(1+delta) and adding
the Chetty mass row.  Bunching mass is the reduced-form observed excess,
governed by {cmd:poolmass} / {cmd:splitmass}.

{pmore}
{cmd:4} -- Saez three-region trapezoid.  Intrinsically a two-point method,
so its bias against a curved h0 comes from evaluating the exact h0 over the
reference regions.

{phang}
{cmd:zstar(}{it:#}{cmd:)} is the kink point z*; all endpoints are centred on it.

{phang}
{cmd:t0(}{it:#}{cmd:)} / {cmd:t1(}{it:#}{cmd:)} are the lower / higher tax
rates; {it:tau = (1 - t0)/(1 - t1)}.

{phang}
{cmd:relslope(}{it:#}{cmd:)} -- degree-1 counterfactual: h0(z*) = 1 and
relative slope h0'(z*)/h0(z*) = {it:#}.  Mutually exclusive with
{cmd:h0poly()}.

{phang}
{cmd:h0poly(}{it:numlist}{cmd:)} -- the counterfactual as a polynomial in
{it:s} = z - z* (levels) or ln z - ln z* (logs), {it:lowest order first}:
{cmd:h0poly(}{it:b0 b1 b2}{cmd:)} means h0(s) = b0 + b1 s + b2 s^2.  Any
degree.  Mutually exclusive with {cmd:relslope()}.

{phang}
{cmd:elasticity(}{it:#}{cmd:)} is the true elasticity against which bias is
evaluated.

{phang}
{cmd:zlo(}{it:#}{cmd:)} / {cmd:zhi(}{it:#}{cmd:)} are the lower / upper bounds
of the estimation (polynomial-fitting) window.  Post-estimation mode reads
the compressed observed bounds from {cmd:e(zlo)} / {cmd:e(zhi)}.

{phang}
{cmd:zl(}{it:#}{cmd:)} / {cmd:zh(}{it:#}{cmd:)} are the lower / upper edges of
the excluded (bunching) region.

{phang}
{cmd:log} calculates in log z.  Post-estimation mode reads this from
{cmd:e(log)}.  As in {helpb polbunch}, the running variable, cutoff, window
and h0 coefficients are assumed to be in ln-earnings already.

{phang}
{cmd:constant} / {cmd:exact} choose how bunching mass is inverted to the
marginal response: {cmd:constant} divides the mass by the fitted density at
z*; {cmd:exact} solves int_0^r h0hat(s) ds = Bhat*bw.  Default: {cmd:constant}
for estimators 1 and 2, {cmd:exact} for 0 and 3.  Match this to how
{cmd:polbunch} was run.

{phang}
{cmd:poolmass} / {cmd:splitmass} choose how bunching mass is backed out:
{cmd:poolmass} subtracts the fitted counterfactual over the whole excluded
window (Chetty 2011 eq. 16); {cmd:splitmass} splits at the kink and uses the
estimator's own h1 above it.  Default: {cmd:poolmass} for estimators 1 and 2,
{cmd:splitmass} for 0 and 3 (estimator 1 is always {cmd:poolmass}).

{phang}
{cmd:bw(}{it:#}{cmd:)} is the bin width used to scale bunching mass; default
{cmd:e(bw_orig)} post-estimation, 1 standalone.

{phang}
{cmd:iterate} solves the self-consistent fixed point
theta_hat = theta_true + bias(theta_true) for the true elasticity and h0
shape (the fitted values are the starting guess).  Useful in post-estimation
mode, where the input elasticity and h0 are themselves biased; not needed
when the inputs are a known DGP.

{phang}
{cmd:tolerance(}{it:#}{cmd:)}, {cmd:maxiter(}{it:#}{cmd:)},
{cmd:underrelax(}{it:#}{cmd:)} control {cmd:iterate}: convergence tolerance
(default 1e-10), iteration cap (default 100), and the under-relaxation
factor that damps the fixed-point step (default 0.5).


{marker examples}{...}
{title:Examples}

{pstd}
Post-estimation mode (immediately after polbunch):

{phang2}{cmd:. polbunch income, cutoff(50000) bw(1000) poly(3) estimator(1) t0(0.2) t1(0.6)}{p_end}
{phang2}{cmd:. polbunchbias}

{pstd}
Post-estimation, overriding the counterfactual with an explicit cubic:

{phang2}{cmd:. polbunchbias, h0poly(1 -0.3 0.1 0.04)}

{pstd}
Standalone, estimator 1, a curved (degree-3) counterfactual, exact inversion:

{phang2}{cmd:. polbunchbias, estimator(1) zstar(1) t0(0.2) t1(0.6) elasticity(0.4) zlo(0) zhi(2.2) zl(0.9) zh(1.15) h0poly(1 -0.5 -0.75 0.25) exact}

{pstd}
Standalone, degree-1 counterfactual (reproduces the earlier behaviour):

{phang2}{cmd:. polbunchbias, estimator(1) zstar(1) t0(0.2) t1(0.6) elasticity(0.4) zlo(0) zhi(2) zl(0.99) zh(1) relslope(0.5) constant}

{pstd}
Estimator 2 with the split-at-the-kink bunching mass and exact inversion:

{phang2}{cmd:. polbunchbias, estimator(2) zstar(1) t0(0.2) t1(0.6) elasticity(0.4) zlo(0) zhi(2) zl(0.99) zh(1) relslope(0.5) splitmass exact}

{pstd}
Saez-style estimator:

{phang2}{cmd:. polbunchbias, estimator(4) zstar(1) t0(0.2) t1(0.6) elasticity(0.4) zlo(0) zhi(2.2) zl(0.99) zh(1) relslope(0)}

{pstd}
Iterated (self-consistent) bias correction:

{phang2}{cmd:. polbunchbias, estimator(1) zstar(1) t0(0.2) t1(0.6) elasticity(0.4) zlo(0) zhi(2) zl(0.99) zh(1) relslope(0.5) iterate}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:polbunchbias} is an {cmd:rclass} command and stores the following in
{cmd:r()}:

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Scalars}{p_end}
{synopt :{cmd:r(estimator)}}estimator code{p_end}
{synopt :{cmd:r(bmodel)}}vestigial, always 0 (the {cmd:bmodel} option was removed){p_end}
{synopt :{cmd:r(islog)}}1 if {cmd:log} was specified; 0 otherwise{p_end}
{synopt :{cmd:r(zstar)}}bunching point z*{p_end}
{synopt :{cmd:r(t0)}}lower tax rate{p_end}
{synopt :{cmd:r(t1)}}higher tax rate{p_end}
{synopt :{cmd:r(tau)}}tax ratio, (1 - t0)/(1 - t1){p_end}
{synopt :{cmd:r(lambda)}}relative slope z* h0'(z*)/h0(z*) used in the final calculation{p_end}
{synopt :{cmd:r(elasticity)}}elasticity used in the final calculation{p_end}
{synopt :{cmd:r(polynomial)}}degree K of the counterfactual polynomial{p_end}
{synopt :{cmd:r(x)}}gross response factor, tau^elasticity{p_end}
{synopt :{cmd:r(rho)}}log response, log(x){p_end}
{synopt :{cmd:r(Delta)}}level proportional response, x - 1{p_end}
{synopt :{cmd:r(zlo)}}lower support point{p_end}
{synopt :{cmd:r(zhi)}}upper support point{p_end}
{synopt :{cmd:r(zL)}}lower excluded-region endpoint{p_end}
{synopt :{cmd:r(zH)}}upper excluded-region endpoint{p_end}
{synopt :{cmd:r(B)}}true bunching mass under the local density model{p_end}
{synopt :{cmd:r(bias_h)}}bias in the fitted density height at z* (= bias in b0){p_end}
{synopt :{cmd:r(bias_slope)}}bias in the fitted s^1 coefficient (= bias in b1){p_end}
{synopt :{cmd:r(bias_lambda)}}bias in the fitted relative slope z* h0'(z*)/h0(z*){p_end}
{synopt :{cmd:r(bias_B)}}bias in the estimated bunching mass{p_end}
{synopt :{cmd:r(bias_response)}}bias in the marginal response (level or log units per {cmd:log}){p_end}
{synopt :{cmd:r(bias_shift)}}bias in the proportional level shift; missing in logs{p_end}
{synopt :{cmd:r(bias_elasticity)}}bias in the elasticity estimate{p_end}
{synopt :{cmd:r(constant)}}1 if the constant-density inversion was used{p_end}
{synopt :{cmd:r(nosplit)}}1 if {cmd:poolmass} was used{p_end}
{synopt :{cmd:r(input_elasticity)}}elasticity supplied to the command{p_end}
{synopt :{cmd:r(corrected_elasticity)}}elasticity after {cmd:iterate}, else the input value{p_end}
{synopt :{cmd:r(iterations)}}number of iterations performed{p_end}
{synopt :{cmd:r(converged)}}1 if {cmd:iterate} converged; 0 if not; missing otherwise{p_end}
{synopt :{cmd:r(iterate_diverged)}}1 if the {cmd:iterate} self-consistency loop failed to converge and the un-iterated bias was reported instead; 0 otherwise{p_end}

{pstd}
For estimators 1 and 2 a missing {cmd:r(bias_h)} means the fitting window was
too short to identify a degree-K polynomial; lower the order or widen the
window.  For estimator 2 it can also mean the right fitting window
[{it:zh},{it:zhi}] on its own is too short to identify the Chetty excess-mass
ratio {it:delta} (widen it).

{p2col 5 28 32 2: Matrices}{p_end}
{synopt :{cmd:r(b)}}1 x (K+6) row vector of displayed biases: {cmd:b0 b1 ... bK} then {cmd:relative_slope number_bunchers marginal_response shift elasticity}{p_end}
{synopt :{cmd:r(bias_beta)}}1 x (K+1) vector of per-coefficient biases, {cmd:b0 ... bK}{p_end}
{synopt :{cmd:r(corrected_beta)}}1 x (K+1) counterfactual polynomial after {cmd:iterate} (renormalised to b0 = 1); the input polynomial otherwise{p_end}

{pstd}
{cmd:b0}, ..., {cmd:bK} are the biases in the coefficients of h0 written as a
polynomial in the centred running variable {it:s}, with h0(z*) normalised to
1.  {cmd:b0} and {cmd:relative_slope} are scale-free; {cmd:b1}, ..., {cmd:bK}
carry units of 1/z^k in levels.


{title:Remarks}

{pstd}
For a true counterfactual that is a degree-K polynomial in the running
variable, the reported biases are exact.  Residual approximation comes only
from binning (the estimator fits bin counts, the formulas integrate
continuously) and, for a real dataset, from the true h0 not being exactly
polynomial.  The h0 fit itself is unweighted least squares, matching
{cmd:polbunch}.

{pstd}
The {cmd:constant}/{cmd:exact} and {cmd:poolmass}/{cmd:splitmass} choices
must match how {cmd:polbunch} was run, or the predicted and realised biases
will differ by the (large) difference between those two inversions.  In
post-estimation mode the axes are taken from {cmd:e()} when {cmd:polbunch}
records them and from the estimator default otherwise.

{pstd}
Because the command is {cmd:rclass}, results are overwritten by the next
{cmd:rclass} command.  Copy {cmd:r(b)} / {cmd:r(bias_beta)} or the scalars
immediately if they are needed later.


{marker references}{...}
{title:References}

{phang}
Andresen, Martin E. (2026). "A better polynomial bunching estimator", working paper.

{phang}
Kleven, Henrik Jacobsen (2016). "Bunching", {it:Annual Review of Economics}.

{phang}
Saez, Emmanuel (2010). "Do Taxpayers Bunch at Kink Points?", {it:American Economic Journal: Economic Policy}.

{phang}
Chetty, Raj, John N. Friedman, Tore Olsen, and Luigi Pistaferri (2011). "Adjustment Costs, Firm Responses, and Micro vs. Macro Labor Supply Elasticities: Evidence from Danish Tax Records", {it:Quarterly Journal of Economics}.


{title:Suggested citation}

{pstd}
Andresen, Martin E. (2026). "POLBUNCH: Stata module for the polynomial bunching estimator." This version VERSION_DATE.{p_end}

{pstd}
Check your installed version date with:{p_end}

{phang2}{cmd:. which polbunch}{p_end}


{marker author}{...}
{title:Author}

{pstd}Martin Eckhoff Andresen{p_end}
{pstd}University of Oslo{p_end}
{pstd}Department of Economics{p_end}
{pstd}Oslo, Norway{p_end}
{pstd}martin.eckhoff.andresen@gmail.com{p_end}


{marker also_see}{...}
{title:Also see}

{p 4 14 2}
Development version: net install polbunch, from("https://raw.githubusercontent.com/martin-andresen/polbunch/master/"){p_end}

{p 7 14 2}
Help: {helpb polbunchplot}, {helpb polbunchgendata}, {helpb polbunchsim}{p_end}
