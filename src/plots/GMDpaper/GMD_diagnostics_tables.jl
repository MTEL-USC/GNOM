
# Compute the diagnostics only once, unless rediagnose=true
(!isdefined(Main, :fDNd_wtags) || rediagnose) && include("GMD_diagnostics_setup.jl")



# Table 3 in GMD paper (source magnitudes, contributions, and bulk ages)
let
    # Dict for long names
    ktostr = Dict(:dust => "Mineral dust",
                  :volc => "Volcanic ash",
                  :sed => "Sedimentary flux",
                  :river => "Riverine discharge",
                  :gw => "Groundwater discharge",
                  :hydro => "Hydrothermal vents",
                  :tot => "Total")

    # units
    uout = Mmol/yr
    uin = mol/m^3/s

    # make a table
    df = DataFrame(
        Symbol("Symbol")=>[symbol2latex(Symbol("σ_$k")) for k in keys(sources)],
        Symbol("Source type")=>[ktostr[k] for k in keys(sources)],
        Symbol("($uout)")=>[∫dV(sₖ*uin,grd)|>uout|>ustrip for sₖ in collect(sources)],
        Symbol("(%source)")=>[∫dV(sₖ,grd)/∫dV(sources.tot,grd)|>per100|>ustrip for sₖ in collect(sources)],
        Symbol("(pM)")=>[totalaverage(DNdₖ*upreferred(uDNd).|>uDNd,grd)|>ustrip for DNdₖ in collect(DNdₖs)],
        Symbol("(%Nd)")=>[∫dV(DNdₖ,grd)/∫dV(DNd,grd)|>per100|>ustrip for DNdₖ in collect(DNdₖs)],
        Symbol("Γ (yr)")=>[∫dV(DNdₖ,grd)/∫dV(sₖ,grd)*s|>yr|>ustrip for (DNdₖ,sₖ) in zip(collect(DNdₖs), collect(sources))]
    )
    @show df

    # Make a LaTeX table
    formatters = (v,i,j) -> (j ≥ 3) ? string("\$", (v ≥ 10 ? Int : identity)(parse(Float64, sprintf1("%.2g", v))), "\$") : v
    @show pretty_table(df, tf=tf_latex_simple, formatters=formatters, nosubheader=true)
end
