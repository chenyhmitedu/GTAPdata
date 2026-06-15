#==
The purpose of the GTAPdata.jl file is to define the module's code and functions. 
It should not contain commands that manage the package environment itself.
==#

module GTAPdata

import JLD2, CSVtoDIC

function io(input::String, output::String)  

    # Convert CSV files to Dictionaries or vectors
    results = CSVtoDIC.source(input)
    d = results[1]              # Note that the element of d has not expanded its keys to the "fullspace". Using d directly may cause missing key errors.
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
    evd     = CSVtoDIC.fullspace(d["evd"], s["set_i"], s["set_g"], s["set_r"])
    evi     = CSVtoDIC.fullspace(d["evi"], s["set_i"], s["set_g"], s["set_r"])

    set_i   = s["set_i"]
    set_g   = s["set_g"]
    set_r   = s["set_r"]
    set_f   = s["set_f"]
    set_cgi = setdiff(set_g, set_i)
    set_sf  = [:lnd, :fix]
    set_mf  = setdiff(set_f, set_sf)
    set_fe  = [:coa, :gas, :p_c]                                                                   

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

    evom  = Dict((f, r) => sum(vfm[f, j, r] for j ∈ set_i)
                        for f ∈ set_f, r ∈ set_r
    )

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

    #=
    vb0 = Dict(
        r => 
            sum((vxmd[(i, s, r)]*(1-rtxs0[(i, s, r)]) + sum(vtwr[(j, i, s, r)] for j ∈ set_i)) for i ∈ set_i, s ∈ set_r)
            -(sum(vxmd[(i, r, s)] for i ∈ set_i, s ∈ set_r)+sum(vst[(i, r)] for i ∈ set_i))            
        for r ∈ [:BGD]
    )
            =#

    #=
    vb_test = Dict(
    r =>
        sum(
            (vxmd[(i, s, r)] * (1 - rtxs0[(i, s, r)]) +
             sum(vtwr[(j, i, s, r)] for j ∈ set_i))
            for i ∈ set_i, s ∈ set_r if s != r
        )
        -
        (
            # Issue: 
            # When r = :BGD, s = :ROW, this is :BGD's exports to ROW, and vst[(i, r)] is BGD's exports to the global pool.
            # When r = :ROW, s = :BGD, this is :ROW's exports to BGD only (small). But vst[(i, r)] is ROW's exports to the global pool (huge).
            # Therefore, while the 1st combination above is BGD's exports to ROW, the 2nd combination can't be ROW's exports to BGD.
            sum(vxmd[(i, r, s)] for i ∈ set_i, s ∈ set_r if s != r)
            + sum(vst[(i, r)] for i ∈ set_i)
        )
    for r ∈ set_r
    )
    =#    

    vafm = Dict(
        (i, g, r) => vdfm[(i, g, r)]*(1+rtfd0[(i, g, r)])+vifm[(i, g, r)]*(1+rtfi0[(i, g, r)])
        for i ∈ set_i, g ∈ set_g, r ∈ set_r
    )

    vafm0 = Dict(
        (i, g, r) => vdfm[(i, g, r)]+vifm[(i, g, r)]
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

    evfm        = Dict((f, r) => sum(vfm[f, j, r] for j ∈ set_i)
                        for f ∈ set_f, r ∈ set_r
    )

    # eind(i,g,r)			= (evd(i,g,r)+evi(i,g,r))/23.88;
    eind        = Dict((i, g, r) => (evd[(i, g, r)]+evi[(i, g, r)])/23.88
                        for i ∈ set_fe, g ∈ set_g, r ∈ set_r
    )

    # Emissions coefficient (bn-tCO2/EJ)
    epslon      = Dict(
    :coa => 0.1*44/12*0.24686,                  
    :p_c => 0.1*44/12*0.199  ,                  
    :gas => 0.1*44/12*0.137                     
    )
    
    # Declare CSAVE elasticities not in GTAP9data package

    esubn       = Dict(i => 0.5 for i ∈ set_g)
    esubve      = Dict(i => 0.4 for i ∈ set_i)
    esubef      = Dict(i => 1.5 for i ∈ set_g)
    esubf       = Dict(i => 1.0 for i ∈ set_g)
    esubc       = 1.0
    esubi       = 1.0

    # Parameters from the CSAVE with the DMR-Devaragan setting

    
    vxmde = Dict((i, r) => sum(vxmd[i, r, s] for s ∈ set_r)
                        for i ∈ set_i, r ∈ set_r
    )

    vxmdm = Dict((i, r) => sum(vxmd[i, s, r] for s ∈ set_r)
                        for i ∈ set_i, r ∈ set_r
    )

    rtxse = Dict(
        (i, r) => begin
            denom = vxmde[i, r] + vst[i, r]
            denom == 0 ? 0.0 :
                sum(rtxs0[i, r, s] for s ∈ set_r) * vxmde[i, r] / denom
        end
        for i ∈ set_i, r ∈ set_r
    )

    e0  = Dict(
        (i, r) => vxmde[i, r]+vst[i, r]
        for i ∈ set_i, r ∈ set_r
    )

    # rtxsm(i)	= sum(nmr, sum(row, rtxs(i,row,nmr)));	
    rtxsm = Dict(
        (i, r) => sum(rtxs0[i, s, r] for s ∈ set_r)
        for i ∈ set_i, r ∈ set_r
    )

    # vtwrm(j,i)	= sum(nmr, sum(row, vtwr(j,i,row,nmr)));
    vtwrm = Dict(
        (j, i, r) => sum(vtwr[j, i, s, r] for s ∈ set_r)
        for j ∈ set_i, i ∈ set_i, r ∈ set_r
    )

    # m0(i)		= vxmdm(i)*(1-rtxsm(i))+sum(j, vtwrm(j,i));
    m0    = Dict(
        (i, r)  => vxmdm[i, r]*(1-rtxsm[i, r]) + sum(vtwrm[j, i, r] for j ∈ set_i)
        for i ∈ set_i, r ∈ set_r   
    )

    # rtmsm(i)	= sum(nmr, sum(row, rtms(i,row,nmr)));
    rtmsm = Dict(
        (i, r) => sum(rtms0[i, s, r] for s ∈ set_r)
        for i ∈ set_i, r ∈ set_r
    )

    # rtfa(i,g,r)$vafm(i,g,r)		= ((vdfm(i,g,r)*(1+rtfd(i,g,r))+vifm(i,g,r)*(1+rtfi(i,g,r)))/vafm(i,g,r))-1;	
    rtfaa = Dict(
        (i, g, r) => vafm0[(i,g,r)] == 0 ? 0.0 :
            (vdfm[(i,g,r)]*(1+rtfd0[(i,g,r)]) +
            vifm[(i,g,r)]*(1+rtfi0[(i,g,r)])) / vafm0[(i,g,r)] - 1
        for i ∈ set_i, g ∈ set_g, r ∈ set_r
    )

    # vafms(i)	= sum(g, vafm0(i,g));
    vafms = Dict(
        (i, r) => sum(vafm0[(i, g, r)] for g ∈ set_g)
        for i ∈ set_i, r ∈ set_r
    )

    # vdfms(i)	= sum(g, vdfm0(i,g));
    vdfms = Dict(
        (i, r) => sum(vdfm[(i, g, r)] for g ∈ set_g)
        for i ∈ set_i, r ∈ set_r
    )

    # vifms(i)	= sum(g, vifm0(i,g));
    vifms = Dict(
        (i, r) => sum(vifm[(i, g, r)] for g ∈ set_g)
        for i ∈ set_i, r ∈ set_r
    )

    # vafmi		= sum(i, vafm0(i,"i")*(1+rtfa0(i,"i")));
    vafmi = Dict(
        r => sum(vafm0[(i, g, r)]*(1+rtfaa[i, g, r]) for i ∈ set_i, g ∈ [:i])
        for r ∈ set_r
    )

    JLD2.@save  output vdfm vxmd vst rtms0 rtxs0 vifm rtfd0 rtfi0 rto0 vfm rtf0 vtwr esubd esubm esubva set_i set_g set_r set_f set_sf set_mf set_cgi set_fe vdm vom pvxmd pvtwr vtw vim vb vafm vafm0 d etadx esub vxm vhm esubn esubi esubve esubef esubf esubc evfm evom rtxse e0 m0 rtmsm rtfaa vafms vdfms vifms vafmi eind epslon





    return
    
end

end # module GTAPdata
