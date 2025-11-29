#using Pkg
Pkg.develop(url = "https://github.com/chenyhmitedu/CSVtoDIC")
#using CSV, DataFrames, CSVtoDIC

module GTAPdata

import JLD2, CSVtoDIC

function io(input::String, output::String)  

    # Convert CSV files to Dictionaries or vectors
    results = CSVtoDIC.source(input)
    d = results[1]
    s = results[2]

    # Expand to the full set space and use GTAP notation
    vdfm    = CSVtoDIC.fullspace(d["vdfm"], s["set_i"], s["set_g"], s["set_r"])                      
    vxmd    = CSVtoDIC.fullspace(d["vxmd"], s["set_i"], s["set_r"], s["set_r"])
    vst     = CSVtoDIC.fullspace(d["vst"], s["set_i"], s["set_r"])
    rtms0   = CSVtoDIC.fullspace(d["rtms"], s["set_i"], s["set_r"], s["set_r"])
    rtxs0   = CSVtoDIC.fullspace(d["rtxs"], s["set_i"], s["set_r"], s["set_r"])
    vifm    = CSVtoDIC.fullspace(d["vifm"], s["set_i"], s["set_g"], s["set_r"])
    rtfd0   = CSVtoDIC.fullspace(d["rtfd"], s["set_i"], s["set_g"], s["set_r"])
    rtfi0   = CSVtoDIC.fullspace(d["rtfi"], s["set_i"], s["set_g"], s["set_r"])
    rto0    = CSVtoDIC.fullspace(d["rto"], s["set_g"], s["set_r"])
    vfm     = CSVtoDIC.fullspace(d["vfm"], s["set_f"], s["set_g"], s["set_r"])
    rtf0    = CSVtoDIC.fullspace(d["rtf"], s["set_f"], s["set_g"], s["set_r"])
    vtwr    = CSVtoDIC.fullspace(d["vtwr"], s["set_i"], s["set_i"], s["set_r"], s["set_r"])
    esubd   = CSVtoDIC.fullspace(d["esubd"], s["set_i"])
    esubm   = CSVtoDIC.fullspace(d["esubm"], s["set_i"])
    esubva  = CSVtoDIC.fullspace(d["esubva"], s["set_g"])
    set_i   = s["set_i"]
    set_g   = s["set_g"]
    set_r   = s["set_r"]
    set_f   = s["set_f"]
    set_cgi = setdiff(set_g, set_i)
    set_sf  = [:lnd, :fix]
    set_mf  = setdiff(set_f, set_sf)
    #set_tr  = [:tran]                                                                   

    # Assignment done in GTAPinGAMS
    d["esub"]       = Dict(i => 0 for i ∈ s["set_g"])       # Top-level elasticity of substitution
    d["esubdm"]     = d["esubd"]                                   

    vdm = Dict(
        (i, r) => sum(vdfm[(i, g, r)] for g in set_g)
        for i in set_i, r in set_r
    )

    vom = Dict(
        (i, r) => sum(vxmd[(i, r, s)] for s in set_r) + vdm[(i, r)] + vst[(i, r)]
        for i ∈ set_i, r ∈ set_r
    )

    vom_a = Dict(
        (i, r) => sum((vdfm[(j, i, r)]*(1+rtfd0[(j, i, r)]) + vifm[(j, i, r)]*(1+rtfi0[(j, i, r)]))/(1-rto0[(i, r)]) for j ∈ set_i) 
        for i ∈ set_cgi, r ∈ set_r
    )

    vom = merge(vom, vom_a)

    vdm_a = Dict(
        (i, r) => vom[(i, r)] for i ∈ [:c, :g], r ∈ set_r
    )

    vdm = merge(vdm, vdm_a)

    #test = sum(vom[(i, r)] for i ∈ set_g, r ∈ set_r)
    #key = [k for (k, v) in vom if v == 0]
    #key = [k for (k, v) in vdm if v == 0]
    #key = [k for (k, v) in vtw if v == 0]
    #key = [k for (k, v) in vim if v == 0]

    pvxmd = Dict(
        (i, s, r) => (1 + rtms0[(i, s, r)])*(1 - rtxs0[(i, s, r)]) for i ∈ set_i, s ∈ set_r, r ∈ set_r
    )

    pvtwr = Dict(
        (i, s, r) => 1 + rtms0[(i, s, r)] for i ∈ set_i, s ∈ set_r, r ∈ set_r
    )

    vtw = Dict(
        j => sum(vst[(j, r)] for r ∈ set_r)
        for j ∈ set_i
    )

    vim = Dict(
        (i, r) => sum(vifm[(i, g, r)] for g ∈ set_g)
        for i ∈ set_i, r ∈ set_r
    )

    vb = Dict(
        r => sum(vom[(i, r)] for i ∈ set_cgi)
            -sum(vfm[(f, g, r)] for f ∈ set_f, g ∈ set_g)
            -sum(vom[(g, r)]*rto0[(g, r)] for g ∈ set_g)
            -sum(vdfm[(i, g, r)]*rtfd0[(i, g, r)] + vifm[(i, g, r)]*rtfi0[(i, g, r)] for i ∈ set_i, g ∈ set_g)
            -sum(vfm[(f, g, r)]*rtf0[(f, g, r)] for f ∈ set_f, g ∈ set_g)
            -sum(rtms0[(i, s, r)]*(vxmd[(i, s, r)]*(1-rtxs0[(i, s, r)]) + sum(vtwr[(j, i, s, r)] for j ∈ set_i)) for i ∈ set_i, s ∈ set_r)
            +sum(rtxs0[(i, r, s)]*vxmd[(i, r, s)] for i ∈ set_i, s ∈ set_r)
        for r ∈ set_r
    )

    vafm = Dict(
        (i, g, r) => vdfm[(i, g, r)]*(1+rtfd0[(i, g, r)])+vifm[(i, g, r)]*(1+rtfi0[(i, g, r)])
        for i ∈ set_i, g ∈ set_g, r ∈ set_r
    )

    etadx = Dict(
        i => esubd[i]*0.5
        for i ∈ set_i
    )

    esub = Dict(
        i => 0.5
        for i ∈ set_i
    )

    merge!(esub, Dict(i => 0.0 for i ∈ set_cgi))

    #vxm(i,s)	= (vst(i,s) + sum(r,vxmd(i,s,r)))
    vxm = Dict((i, s) => vst[i, s] + sum(vxmd[i, s, r] for r ∈ set_r)
               for i ∈ set_i, s ∈ set_r 
    )

    vhm = Dict((j, r) => vom[j, r]-vxm[j, r]
               for j ∈ set_i, r ∈ set_r 
    )

    for i ∈ set_cgi, r ∈ set_r 
        vom[i, r] = sum((vdfm[j, i, r]*(1 + rtfd0[j, i, r]) + vifm[j, i, r]*(1 + rtfi0[j, i, r]))/(1 - rto0[i, r]) for j ∈ set_i)
    end

    evom        = Dict((f, r) => sum(vfm[f, j, r] for j ∈ set_i)
                        for f ∈ set_f, r ∈ set_r
    )

    # Declare CSAVE elasticities not in GTAP9data package

    esubn       = Dict(i => 0.5 for i ∈ set_g)
    esubve      = Dict(i => 0.4 for i ∈ set_i)
    esubef      = Dict(i => 1.5 for i ∈ set_g)
    esubf       = Dict(i => 1.0 for i ∈ set_g)
    esubc       = 1.0

    JLD2.@save output vdfm vxmd vst rtms0 rtxs0 vifm rtfd0 rtfi0 rto0 vfm rtf0 vtwr esubd esubm esubva set_i set_g set_r set_f set_sf set_mf set_cgi vdm vom pvxmd pvtwr vtw vim vb d vafm etadx esub vxm vhm esubn esubve esubef esubf esubc evom

    return IO
end

end # module GTAPdata
