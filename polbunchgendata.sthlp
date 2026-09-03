{smcl}
{cmd:help polbunchgendata}
{hline}

{title:Title}

{p2colset 5 14 16 2}{...}
{p2col:{cmd:polbunchgendata} {hline 2}}Simulates data from the iso-elastic labor supply model in a setting with a kink.{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 15 2}
{cmd:polbunchgendata} {it:newvar} {cmd:,} [{it:options}]

{pstd}
{cmd:polbunchgendata} generates individual-level earnings from the iso-elastic
labor supply model with a convex tax kink at {opt cutoff()} (marginal rate
{opt t0()} below, {opt t1()} above).  {opt distribution()} gives the
{it:counterfactual} earnings density -- what people would earn facing {opt t0()}
everywhere.  Bunchers are moved to the kink; non-bunchers above the kink are
relocated according to the model.  The behavioural response is written into
{it:newvar}.  The current dataset is replaced.

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt obs(#)}}Number of observations.  Default 5000.{p_end}
{synopt:{opt cutoff(#)}}Kink point (in levels, or in logs with {opt log}).  Default 1.{p_end}
{synopt:{opt t0(#)}}Marginal tax rate below the kink.  Default 0.2.{p_end}
{synopt:{opt t1(#)}}Marginal tax rate above the kink.  Default 0.6.{p_end}
{synopt:{opt el(spec)}}Compensated elasticity.  {it:spec} is {bf:either} a
nonnegative number (homogeneous elasticity) {bf:or} any Stata expression,
evaluated per observation, that draws an individual elasticity -- e.g.
{cmd:el(0.4)}, {cmd:el(0.4 + rnormal(0,0.1))}, {cmd:el(rbeta(2,3)*1.5)}.
Draws are floored at 1e-8.  Default 0.4.{p_end}
{synopt:{opt incomeeffect(spec)}}Income effects.  {it:spec} is {bf:either} a
number in [0,1) {bf:or} a Stata expression, evaluated per observation, drawing
an individual curvature parameter eta_i (values outside [0,1) are clipped, with
a note).  eta is the consumption-curvature parameter of
u(c,z) = c^(1-eta)/(1-eta) - (n/(1+1/e))(z/n)^(1+1/e); the h(z) curvature
exponent is 1/e_i = 1/ec_i - eta_i.  eta=0 (the default) reproduces the
no-income-effect iso-elastic model exactly.  The bunching window is unchanged
(Saez 2010); only the relocation of non-bunchers above the kink differs, and it
is solved individually.  Not allowed with {opt log}.{p_end}
{synopt:{opt distribution(string)}}Expression generating the counterfactual
earnings distribution, e.g. {cmd:rbeta(2,5)}.  {cmd:triangular(a,b,c)} is also
allowed, with {cmd:triangular(0,3,0)} the default.{p_end}
{synopt:{opt log}}Earnings (and {opt cutoff()}) are in logs.{p_end}
{synopt:{opt buncherror(op expr)}}Optimisation friction / measurement error for
bunchers, applied as {cmd:cutoff}{it: op expr}.  Give an operator followed by an
expression, e.g. {cmd:buncherror(+rnormal(0,0.05))} or
{cmd:buncherror(*exp(rnormal(0,0.1)))}.{p_end}
{synoptline}


{marker examples}{...}
{title:Examples}

{pstd}Homogeneous elasticity, default triangular counterfactual:{p_end}
{phang2}{cmd:. polbunchgendata z, obs(20000) cutoff(1) el(0.4) t0(0.2) t1(0.6)}{p_end}

{pstd}Heterogeneous elasticity centred on 0.4:{p_end}
{phang2}{cmd:. polbunchgendata z, el(0.4 + rnormal(0,0.15)) t0(0.2) t1(0.6)}{p_end}

{pstd}Heterogeneous income effects and an optimisation friction:{p_end}
{phang2}{cmd:. polbunchgendata z, el(0.5) incomeeffect(rbeta(2,3)*0.6) buncherror(+rnormal(0,0.03))}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}{cmd:polbunchgendata} stores the following in {cmd:r()}:{p_end}
{synoptset 20 tabbed}{...}
{synopt:{cmd:r(n_bunchers)}}number of individuals who bunch at the kink{p_end}
{synopt:{cmd:r(share_bunching)}}{cmd:r(n_bunchers)} / {opt obs()}{p_end}
{synopt:{cmd:r(el_mean)}}realised mean of the drawn compensated elasticity{p_end}
{synopt:{cmd:r(incomeeffect)}}realised mean of the drawn income-effect parameter{p_end}
{synoptline}


{marker Author}{...}
{title:Author}

{pstd}Martin Eckhoff Andresen{p_end}
{pstd}University of Oslo{p_end}
{pstd}Department of Economics{p_end}
{pstd}Oslo, Norway{p_end}
{pstd}martin.eckhoff.andresen@gmail.com{p_end}

{marker also_see}{...}
{title:Also see}

{p 7 14 2}
{helpb polbunch} for help on the main command {cmd:polbunch}.{p_end}
