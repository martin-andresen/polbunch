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
{cmd:polbunchplot} [namelist] [{cmd:,} {opt names(string)} {opt graph_opts(string)} {opt leg:end_opts(string)} {opt limit(numlist)} {opt log} {opt tru:ncate}]


{pstd}
{cmd:polbunchplot} plots bunching plots after polbunch estimation for the running variable, based on the polbunch estimate stored in [namelist], if specified, or in memory.


{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt names(string)}} Custom legend labels for the plotted models, one per model in the same order as {it:namelist}, separated by a vertical bar {cmd:|}. Labels may contain spaces and punctuation. Applies to the multiple-model legend only; the number of labels must match the number of models. Example: {cmd:names("Chetty et al. (2011)|Naive polynomial|This paper")}.{p_end}
{synopt:{opt graph_opts(string)}} Options passed through to the final {helpb twoway} call (e.g. {cmd:name()}, {cmd:title()}, {cmd:xtitle()}, {cmd:scheme()}).{p_end}
{synopt:{opt leg:end_opts(string)}} Suboptions passed through to {helpb legend_option:legend()} (e.g. {cmd:cols()}, {cmd:pos()}, {cmd:region()}, {cmd:label()}). To rename model keys, prefer {opt names()}: legend keys are numbered by plot, not by model, so {cmd:label(#)} is hard to target.{p_end}
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
