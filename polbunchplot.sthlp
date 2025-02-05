{smcl}
{cmd:help polbunchplot}
{hline}

{title:Title}

{p2colset 5 14 16 2}{...}
{p2col:{cmd:polbunchplot} {hline 2}}Plotting of bunching estimates from {cmd: polbunch}
{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 10 15 2}
{cmd:polbunchplot} [namelist] [{cmd:,} graph_opts(string) limit(numlist) log truncate]


{pstd}
{cmd:polbunchplot} plots bunching plots after polbunch estimation for the running variable, based on the polbunch estimate stored in [namelist], if specified, or in memory.


{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt limit(numlist)}} Only plot values of earnings between the two numbers in limit().{p_end}
{synopt:{opt tru:ncate}} Truncate values in the bunching region to be no larger than the maximum outside of the bunching region; useful if bunching is so substantial the figure cannot be used to evaluate fit.{p_end}
{synopt:{opt log}} present log frequency on the y axis.{p_end}
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
{helpb polbunch} for help on the main command {cmd: polbunch}.{p_end}
