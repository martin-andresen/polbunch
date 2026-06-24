{smcl}
{* *! version 1.0.0 23jun2026}{...}
{vieweralsosee "[R] return" "help return"}{...}
{viewerjumpto "Syntax" "polbunchbias##syntax"}{...}
{viewerjumpto "Description" "polbunchbias##description"}{...}
{viewerjumpto "Options" "polbunchbias##options"}{...}
{viewerjumpto "Examples" "polbunchbias##examples"}{...}
{viewerjumpto "Stored results" "polbunchbias##results"}{...}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{hi:polbunchbias} {hline 2}}Calculate first-order bias of polynomial bunching estimators{p_end}
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
[{cmd:lambda(}{it:#}{cmd:)}]
[{it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required (standalone mode only)}
{synopt :{cmd:estimator(}{it:#}{cmd:)}}estimator whose first-order bias is evaluated{p_end}
{synopt :{cmd:zstar(}{it:#}{cmd:)}}bunching point or kink point, z*{p_end}
{synopt :{cmd:t0(}{it:#}{cmd:)}}lower tax rate, t0{p_end}
{synopt :{cmd:t1(}{it:#}{cmd:)}}higher tax rate, t1{p_end}
{synopt :{cmd:elasticity(}{it:#}{cmd:)}}elasticity used to generate the implied behavioral response{p_end}
{synopt :{cmd:zlo(}{it:#}{cmd:)}}lower support point used in the bias calculation{p_end}
{synopt :{cmd:zhi(}{it:#}{cmd:)}}upper support point used in the bias calculation{p_end}
{synopt :{cmd:zl(}{it:#}{cmd:)}}lower edge of the excluded or bunching region{p_end}
{synopt :{cmd:zh(}{it:#}{cmd:)}}upper edge of the excluded or bunching region{p_end}

{syntab:Optional (both modes)}
{synopt :{cmd:lambda(}{it:#}{cmd:)}}relative slope of the counterfactual density at z*; in post-estimation mode, computed automatically from {cmd:e(b)} if omitted{p_end}
{synopt :{cmd:bmodel(}{it:#}{cmd:)}}report model-implied bunching response for estimators that support it; default is {cmd:bmodel(0)}{p_end}
{synopt :{cmd:log}}perform the calculation in log z rather than level z (standalone mode only; post-estimation mode reads this from {cmd:e(log)}){p_end}
{synopt :{cmd:bw(}{it:#}{cmd:)}}bin width used to scale bunching mass; default is {cmd:e(bw_orig)} in post-estimation mode and {cmd:bw(1)} in standalone mode{p_end}
{synopt :{cmd:iterate}}iteratively correct the input elasticity and slope by subtracting the estimated bias{p_end}
{synopt :{cmd:tolerance(}{it:#}{cmd:)}}convergence tolerance for {cmd:iterate}; default is {cmd:tolerance(1e-10)}{p_end}
{synopt :{cmd:maxiter(}{it:#}{cmd:)}}maximum number of iterations for {cmd:iterate}; default is {cmd:maxiter(100)}{p_end}
{synopt :{cmd:constant}}use the constant-density approximation when translating bunching mass into the marginal response{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:polbunchbias} calculates analytical first-order bias expressions for
polynomial bunching estimators under a local linear counterfactual density.
The command is intended as a diagnostic and approximation tool: it takes the
structural primitives supplied by the user, computes the implied bunching
response, and reports the first-order bias in the estimated density height,
slope, bunching mass, marginal response, and elasticity.

{pstd}
{cmd:polbunchbias} operates in two modes.  When called without any of
the nine required options, it enters {it:post-estimation mode} and reads all
parameters from the most recent {cmd:polbunch} results stored in {cmd:e()}.
The kink point, bin width, excluded-region endpoints, estimation window,
log/level flag, tax rates, and elasticity are all read automatically.
For polynomial estimators (1–3), the relative slope lambda is derived from
the first-order coefficient in {cmd:h0}.  For the Saez estimator (4), lambda
is imputed from the implicit two-point counterfactual defined by
{cmd:_b[h0:_cons]} and {cmd:_b[h1:_cons]}.  The user may override lambda by
supplying {cmd:lambda()} explicitly.

{pstd}
When all nine required options are supplied, the command operates in
{it:standalone mode} and uses only the supplied values.  Supplying a subset
of the nine required options is an error.

{pstd}
The command assumes a normalized counterfactual density with height one at
the bunching point. The local relative slope is supplied by {cmd:lambda()}.
Internally, the density slope is set to {it:m=lambda/zstar}. The tax change
implies

{p 12 12 2}
{it:tau} = (1 - {it:t0})/(1 - {it:t1})

{pstd}
and the response is determined by the supplied elasticity. In levels,
the proportional shift is {it:Delta = tau^elasticity - 1}; in logs, the
corresponding log response is {it:rho = log(tau^elasticity)}.

{pstd}
The reported biases are stored in {cmd:r(b)} and are also displayed as a
one-row table with seven columns:
{cmd:h}, {cmd:slope}, {cmd:relative_slope}, {cmd:number_bunchers},
{cmd:marginal_response}, {cmd:shift}, and {cmd:elasticity}.

{pstd}
The current implementation supports estimator codes {cmd:1}, {cmd:2}, and
{cmd:4}. Unsupported estimator codes return missing values.


{marker options}{...}
{title:Options}

{phang}
{cmd:estimator(}{it:#}{cmd:)} specifies the estimator whose first-order bias
is evaluated. The current implementation supports:

{pmore}
{cmd:estimator(1)}: polynomial estimator that estimates the counterfactual
density from the included left and right regions and then computes bunching
mass as actual mass in the excluded region minus the fitted counterfactual
mass.

{pmore}
{cmd:estimator(2)}: Chetty-style estimator imposing a restriction between
pre- and post-reform densities. In levels, the estimator also incorporates
a mass-row restriction. With {cmd:bmodel(0)}, the final reported bunching
mass is reduced-form; with {cmd:bmodel(1)}, the final response is
model-implied.

{pmore}
{cmd:estimator(4)}: Saez-style three-region trapezoid approximation using
left and right reference densities and an implicit transformation from
bunching mass to the behavioral response.

{phang}
{cmd:zstar(}{it:#}{cmd:)} specifies the kink or notch point z*. All support
and excluded-region endpoints are centered around this point internally.

{phang}
{cmd:t0(}{it:#}{cmd:)} and {cmd:t1(}{it:#}{cmd:)} specify the lower and higher
tax rates. The tax ratio is computed as
{it:tau = (1 - t0)/(1 - t1)}.

{phang}
{cmd:lambda(}{it:#}{cmd:)} specifies the relative slope of the counterfactual
density at z*. With the normalization h0(z*) = 1, the level slope is
{it:m = lambda/zstar}.

{phang}
{cmd:elasticity(}{it:#}{cmd:)} specifies the elasticity used to generate the
true response against which first-order bias is evaluated.

{phang}
{cmd:zlo(}{it:#}{cmd:)} and {cmd:zhi(}{it:#}{cmd:)} specify the lower and upper
support points used in the calculation.

{phang}
{cmd:zl(}{it:#}{cmd:)} and {cmd:zh(}{it:#}{cmd:)} specify the lower and upper
edges of the excluded or bunching region.

{phang}
{cmd:bmodel(}{it:#}{cmd:)} controls how estimator 2 reports the final
bunching response. The default, {cmd:bmodel(0)}, reports the reduced-form
bunching-mass calculation. {cmd:bmodel(1)} reports the model-implied
response where supported.

{phang}
{cmd:log} requests the log-z version of the calculation. Without {cmd:log},
the calculation is performed in levels.

{phang}
{cmd:bw(}{it:#}{cmd:)} specifies the bin width used to scale bunching mass.
The default is {cmd:bw(1)}.

{phang}
{cmd:iterate} iteratively subtracts the estimated bias from the input
elasticity and slope parameters and recomputes the bias at the corrected
values. Iteration updates both elasticity and lambda.

{phang}
{cmd:tolerance(}{it:#}{cmd:)} specifies the convergence tolerance for
{cmd:iterate}. The default is {cmd:tolerance(1e-10)}.

{phang}
{cmd:maxiter(}{it:#}{cmd:)} specifies the maximum number of iterations for
{cmd:iterate}. The default is {cmd:maxiter(100)}.

{phang}
{cmd:constant} uses the constant-density approximation when translating
bunching mass into the marginal response. Without {cmd:constant},
{cmd:polbunchbias} solves the quadratic equation implied by the local linear
density.


{marker examples}{...}
{title:Examples}

{pstd}
Post-estimation mode (immediately after polbunch):

{phang2}{cmd:. polbunch income, cutoff(50000) bw(1000) poly(3) estimator(1) t0(0.2) t1(0.6)}{p_end}
{phang2}{cmd:. polbunchbias}

{pstd}
Post-estimation mode, overriding lambda:

{phang2}{cmd:. polbunchbias, lambda(0.3)}

{pstd}
Estimator 1, standalone level calculation:

{pstd}
Estimator 1 using the constant-density approximation:

{phang2}{cmd:. polbunchbias, estimator(1) zstar(1) t0(0.2) t1(0.6) lambda(0.5) elasticity(0.4) zlo(0) zhi(2) zl(0.99) zh(1) constant}

{pstd}
Estimator 2 with the model-implied response:

{phang2}{cmd:. polbunchbias, estimator(2) bmodel(1) zstar(1) t0(0.2) t1(0.6) lambda(0.5) elasticity(0.4) zlo(0) zhi(2) zl(0.99) zh(1)}

{pstd}
Saez-style estimator:

{phang2}{cmd:. polbunchbias, estimator(4) zstar(1) t0(0.2) t1(0.6) lambda(0) elasticity(0.4) zlo(0) zhi(2.2) zl(0.99) zh(1)}

{pstd}
Iterated bias correction:

{phang2}{cmd:. polbunchbias, estimator(1) zstar(1) t0(0.2) t1(0.6) lambda(0.5) elasticity(0.4) zlo(0) zhi(2) zl(0.99) zh(1) iterate}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:polbunchbias} is an {cmd:rclass} command and stores the following in
{cmd:r()}:

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Scalars}{p_end}
{synopt :{cmd:r(estimator)}}estimator code{p_end}
{synopt :{cmd:r(bmodel)}}bmodel option value{p_end}
{synopt :{cmd:r(islog)}}1 if {cmd:log} was specified; 0 otherwise{p_end}
{synopt :{cmd:r(zstar)}}bunching point z*{p_end}
{synopt :{cmd:r(t0)}}lower tax rate{p_end}
{synopt :{cmd:r(t1)}}higher tax rate{p_end}
{synopt :{cmd:r(tau)}}tax ratio, (1 - t0)/(1 - t1){p_end}
{synopt :{cmd:r(lambda)}}lambda used in the final calculation{p_end}
{synopt :{cmd:r(elasticity)}}elasticity used in the final calculation{p_end}
{synopt :{cmd:r(x)}}gross response factor, tau^elasticity{p_end}
{synopt :{cmd:r(rho)}}log response, log(x){p_end}
{synopt :{cmd:r(Delta)}}level proportional response, x - 1{p_end}
{synopt :{cmd:r(zlo)}}lower support point{p_end}
{synopt :{cmd:r(zhi)}}upper support point{p_end}
{synopt :{cmd:r(zL)}}lower excluded-region endpoint{p_end}
{synopt :{cmd:r(zH)}}upper excluded-region endpoint{p_end}
{synopt :{cmd:r(B)}}true bunching mass under the local density model{p_end}
{synopt :{cmd:r(bias_h)}}bias in the estimated density height at z*{p_end}
{synopt :{cmd:r(bias_B)}}bias in estimated bunching mass{p_end}
{synopt :{cmd:r(bias_response)}}bias in the marginal response, in level or log units depending on {cmd:log}{p_end}
{synopt :{cmd:r(bias_shift)}}bias in the proportional level shift; missing in log calculations{p_end}
{synopt :{cmd:r(bias_elasticity)}}bias in the elasticity estimate{p_end}
{synopt :{cmd:r(bias_slope)}}bias in the estimated density slope{p_end}
{synopt :{cmd:r(bias_lambda)}}bias in the estimated relative slope{p_end}
{synopt :{cmd:r(constant)}}1 if {cmd:constant} was specified; 0 otherwise{p_end}
{synopt :{cmd:r(input_elasticity)}}elasticity originally supplied by the user{p_end}
{synopt :{cmd:r(input_lambda)}}lambda originally supplied by the user{p_end}
{synopt :{cmd:r(corrected_elasticity)}}final elasticity after iteration, or input elasticity if not iterated{p_end}
{synopt :{cmd:r(corrected_lambda)}}final lambda after iteration, or input lambda if not iterated{p_end}
{synopt :{cmd:r(iterations)}}number of iterations performed{p_end}
{synopt :{cmd:r(converged)}}1 if iteration converged; 0 if not; missing if {cmd:iterate} was not specified{p_end}

{p2col 5 28 32 2: Matrices}{p_end}
{synopt :{cmd:r(b)}}1 x 7 row vector of displayed biases{p_end}

{pstd}
The columns of {cmd:r(b)} are:

{p2colset 9 30 32 2}{...}
{p2col :{cmd:h}}bias in estimated density height at z*{p_end}
{p2col :{cmd:slope}}bias in estimated density slope{p_end}
{p2col :{cmd:relative_slope}}bias in estimated relative slope lambda{p_end}
{p2col :{cmd:number_bunchers}}bias in estimated bunching mass{p_end}
{p2col :{cmd:marginal_response}}bias in the marginal response{p_end}
{p2col :{cmd:shift}}bias in the proportional level shift{p_end}
{p2col :{cmd:elasticity}}bias in the elasticity estimate{p_end}
{p2colreset}


{title:Remarks}

{pstd}
The command evaluates bias formulas at the user-supplied local linear
counterfactual density. The calculations are first-order approximations to
the behavior of the corresponding bunching estimators. The {cmd:constant}
option applies a constant-density approximation in the final transformation
from bunching mass to the marginal response; otherwise, the command solves
the quadratic equation implied by the local linear density.

{pstd}
Because the command is {cmd:rclass}, results are overwritten by the next
{cmd:rclass} command. Use {cmd:matrix b = r(b)} or save the returned scalars
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
