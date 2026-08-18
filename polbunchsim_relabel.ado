*! polbunchsim_relabel - clean up simulate's auto-generated result names,
*! then (by default) reshape the results to long format.
*!
*! For a single btype()/estimator()/clist() combination (numest==1),
*! polbunchsim keeps polbunch's own equation names with no model prefix,
*! so btype/estimator/cval must be supplied explicitly to include those
*! columns in the reshape: polbunchsim_relabel, btype(1) estimator(1) cval(0)
program define polbunchsim_relabel
    syntax [, NOLong BType(string) ESTimator(string) CVal(string)]

    /*
        simulate collects e(b) columns by default and names each result
        variable itself, capped at 32 chars. When the eq:coef pair doesn't
        fit, simulate falls back to a generic name (_sim_43 etc) -- but it
        still attaches a "[eq]_b[coef]" variable label to every collected
        variable, whether or not the name itself was descriptive. Parse
        that label back rather than guessing simulate's internal 32-char
        budget (it inserts an undocumented literal "_b_" token between eq
        and coef that's easy to get wrong).

        polbunchsim's equation names follow a fixed pattern of its own:
        b<btype>e<estimator>c<cval>[_<subeq>], where <subeq> is polbunch's
        own equation (h0, h1, bunching, ...). That pattern -- not the
        number of underscores in a coefficient's name -- is what drives
        the optional reshape to long below, so coefficient names like
        "number_bunchers" or "marginal_response" don't need special-casing
        the way they would with plain string-splitting.
    */
    quietly ds
    local allvars `r(varlist)'
    local used `allvars'

    local nrenamed = 0
    local k = 0

    foreach v of local allvars {
        local lbl : variable label `v'
        if `"`lbl'"' == "" continue
        if !regexm(`"`lbl'"', "^\[?([^]]*)\]?_b\[(.+)\]$") continue

        local eq   = regexs(1)
        local coef = regexs(2)

        if "`eq'" == "" | "`eq'" == "b" local base "`coef'"
        else local base "`eq'_`coef'"
        local base = subinstr("`base'", " ", "_", .)

        /*
            polbunch's polynomial terms use factor-variable notation
            (c.z#c.z, c.z#c.z#c.z, ...) which is valid in a matrix
            coefficient name but not in a Stata variable name: "." and
            "#" would make -rename- error out. Strip the "c." operator
            prefix and turn "#" into "_" so the name stays readable and
            valid; the untouched original text is still kept as the
            -stat- value if this gets reshaped to long below.
        */
        local base = subinstr("`base'", "c.", "", .)
        local base = subinstr("`base'", "i.", "", .)
        local base = subinstr("`base'", "b.", "", .)
        local base = subinstr("`base'", "o.", "", .)
        local base = subinstr("`base'", "#", "_", .)
        local base = subinstr("`base'", ".", "_", .)

        local newname = substr("`base'", 1, 32)
        if "`newname'" != "`v'" {
            local sfx_k = 1
            while strpos(" `used' ", " `newname' ") > 0 {
                local sfx "_`sfx_k'"
                local room = 32 - strlen("`sfx'")
                local newname = substr("`base'", 1, `room') + "`sfx'"
                local ++sfx_k
            }
            local used = subinstr(" `used' ", " `v' ", " ", .)
            local used `used' `newname'
            rename `v' `newname'
            local ++nrenamed
        }
        label variable `newname' `"`lbl'"'

        local ++k
        local var`k' "`newname'"
        local coef`k' "`coef'"
        local bt`k'  ""
        local est`k' ""
        local cv`k'  ""
        local sub`k' ""
        local hasprefix`k' = 0
        if regexm("`eq'", "^b([0-9]+)e([0-9]+)c([0-9]+)(_(.+))?$") {
            local bt`k'  = regexs(1)
            local est`k' = regexs(2)
            local cv`k'  = regexs(3)
            local sub`k' = regexs(5)
            local hasprefix`k' = 1
        }
        else if "`eq'" != "b" {
            local sub`k' "`eq'"
        }
    }

    display as text "polbunchsim_relabel: renamed `nrenamed' variable(s) using their simulate labels"

    if "`nolong'" != "" exit

    /*
        With only one btype()/estimator()/clist() combination, polbunchsim
        never builds a b<btype>e<estimator>c<cval> prefix at all -- it just
        keeps polbunch's own equation names (h0, h1, bunching, ...), since
        there's nothing to disambiguate between models. Those columns are
        still reshapable, but btype/estimator/cval are constants that only
        the caller knows (they lived in the polbunchsim syntax, not in any
        e()-return that survives past `simulate`) -- pass them here.
    */
    local have_const = ("`btype'" != "" & "`estimator'" != "" & "`cval'" != "")

    local matched ""
    local skipped = 0
    forvalues i = 1/`k' {
        if `hasprefix`i'' {
            local matched `matched' `i'
        }
        else if `have_const' {
            local bt`i'  "`btype'"
            local est`i' "`estimator'"
            local cv`i'  "`cval'"
            local matched `matched' `i'
        }
        else local ++skipped
    }

    if "`matched'" == "" {
        if `skipped' > 0 {
            display as error "polbunchsim_relabel: found `skipped' column(s) with no b<btype>e<estimator>c<cval> prefix (a single-combo simulate run). Re-run with btype()/estimator()/cval() to supply the fixed values, e.g. polbunchsim_relabel, btype(1) estimator(1) cval(0)"
        }
        else {
            display as text "polbunchsim_relabel: no b<btype>e<estimator>c<cval>-style columns found; skipping reshape"
        }
        exit
    }
    if `skipped' > 0 {
        display as text "polbunchsim_relabel: `skipped' column(s) with no model prefix were skipped (pass btype()/estimator()/cval() to include them)"
    }

    capture confirm variable simid
    if _rc gen long simid = _n

    tempfile _long
    local first = 1
    foreach i of local matched {
        preserve
            keep simid `var`i''
            rename `var`i'' value
            gen btype     = `bt`i''
            gen estimator = `est`i''
            gen cval      = `cv`i''
            gen subeq     = "`sub`i''"
            gen stat      = "`coef`i''"
            if `first' {
                save `_long', replace
            }
            else {
                append using `_long'
                save `_long', replace
            }
        restore
        local first = 0
    }

    use `_long', clear
    order simid btype estimator cval subeq stat value
    sort simid btype estimator cval subeq stat

    local nmatched : word count `matched'
    display as text "polbunchsim_relabel: reshaped `nmatched' matched column(s) to long; other columns (if any) were dropped"
end
