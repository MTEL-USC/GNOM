# solve it once
s = solve(prob, CTKAlg(), preprint="Nd & εNd solve ", τstop=ustrip(u"s", 1e3u"Myr"))

#================================================
Optimization setup
================================================#

# modifier function for objective function
# (takes tracers as input and spits out composites
#  that are compared to observations).
modify(DSi, PSi) = (DSi,)

# Weights for objective function
const ωs = (1.0,) # the weight for the mismatches
const ωp = 1e-4   # parameter weights (small)
const Si_obs_tuple = let
    obs = WorldOceanAtlasTools.observations("silicate")
    obs.value = ustrip.(upreferred.(obs.silicate * ρSW))
    (obs,)
end

f, ∇ₓf = f_and_∇ₓf(ωs, ωp, grd, modify, Si_obs_tuple, SiParams)

# Use F1 for gradient and Hessian
λ = p2λ(p)
τstop = ustrip(u"s", 1e3u"Myr")
mem = F1Method.initialize_mem(F, ∇ₓf, ∇ₓF, x0, λ, CTKAlg(); preprint="mem ", τstop=τstop)

function objective(λ)
    p = λ2p(SiParams, λ) ; @show p
    F1Method.objective(f, F, ∇ₓF, mem, λ, CTKAlg(), preprint="obj ", τstop=τstop)
end
gradient(λ) = F1Method.gradient(f, F, ∇ₓf, ∇ₓF, mem, λ, CTKAlg(), preprint="grad", τstop=τstop)
hessian(λ) = F1Method.hessian(f, F, ∇ₓf, ∇ₓF, mem, λ, CTKAlg(), preprint="hess ", τstop=τstop)

# Reduced g_tol for optimization
opt = Optim.Options(store_trace=false, show_trace=true, extended_trace=false, g_tol=1e-3)

#================================================
Optimization run
================================================#
results = optimize(objective, gradient, hessian, λ, NewtonTrustRegion(), opt; inplace=false)

p_optimized = λ2p(SiParams, results.minimizer)
prob_optimized = SteadyStateProblem(fun, x0, p_optimized)
s_optimized = solve(prob_optimized, CTKAlg(), τstop=ustrip(u"s", 1e3u"Myr")).u

# TODO find a more generic approach to save this data... Maybe I can use datadeps?

DSi, PSi = unpack_tracers(s_optimized, grd)
DSi, = modify(DSi, PSi)
tp_opt = AIBECS.table(p_optimized)

# create directories to write in
jldsave(joinpath(output_path, "optimized_Simodel_$circname.jld2"); DSi, PSi, s_optimized, tp_opt, headcommit)

println("Done!")
# print optimized parameters in md table
open(joinpath(output_path, "optimized_Si_parameters_$circname.md"), "w") do f
    pretty_table(f, AIBECS.table(p_optimized)[:,[[1,2,3,4];6:end]], nosubheader=true, tf=tf_markdown)
end


# TODO: Make this table into a LaTeX table and print it out for GMD paper

#=
# and print a LaTeX table for paper
function symbol2latex(s)
    s = latexstring(s)
    s = replace(s, "σ" => "\\sigma")
    s = replace(s, "α" => "\\alpha")
    s = replace(s, "ε" => "\\varepsilon")
    s = replace(s, "_dst" => "")
    s = replace(s, "′" => "'") # replace spaces with small spaces
    s = replace(s, r"_([A-Za-z]*)" => s"_\\mathrm{\1}")
    s = replace(s, r"₀" => "")
end
function latexify(U)
    str = string(U)
    str = replace(str, "⁻" => s"^-") # replace spaces with small spaces
    str = replace(str, "‱" => s"") # replace spaces with small spaces
    str = replace(str, "¹" => "1") # replace spaces with small spaces
    str = replace(str, "²" => "2") # replace spaces with small spaces
    str = replace(str, r"\^-(?<exp>\d+)" => s"$^{-\g<exp>}$") # add brackets around exponents
    str = replace(str, r"\s" => s"\\,") # replace spaces with small spaces
    return str
end
latexeNd(s) = replace(s, "εNd"=>"\\eNd{}")
latexbool(optimized) = optimized ? "\\checkmark" : ""
tp2 = select(tp_opt,
             :Symbol=>ByRow(symbol2latex)=>:Symbol,
             :Value,
             Symbol("Initial value"),
             :Unit=>ByRow(latexify)=>:Unit,
             :Description=>ByRow(latexeNd)=>:Description,
             :Optimizable=>ByRow(latexbool)=>:Optimized)
open(joinpath(archive_path, "optimized_parameters.tex"), "w") do f
    pretty_table(f, tp2, tf=tf_latex_simple, formatters = ft_printf("%.3g", [2, 3]), nosubheader=true)
end
=#

# TODO: Convert the plots below to Makie and make them quality for GMD paper


#=
# Plot optimal model skill
function _locations(xmodel, obs, ux)
    M = interpolationmatrix(grd, obs)
    xmodelatobs = M * xmodel
    xdepthatobs = M * depthvec(grd)
    iobswet = findall(iswet(grd, obs))
    xwetobs = uconvert.(ux, obs.value[iobswet] * upreferred(ux))
    return xmodelatobs, xdepthatobs, iobswet, xwetobs
end
ux = mmol/m^3
DSi_model = uconvert.(ux, DSi * upreferred(ux))
xmodelatobs, xdepthatobs, iobswet, xwetobs = _locations(DSi_model, Si_obs_tuple[1], ux)

boundary = (-10,210)
x, y = ustrip.(xwetobs), ustrip.(xmodelatobs)
bw = (boundary[2]-boundary[1])/256
D = kde((x, y); boundary=(boundary, boundary), bandwidth=(bw,bw))
# calculate cumulative density from density
δx = step(D.x)
δy = step(D.y)
Q = vec(D.density) * δx * δy
idx = sortperm(Q)
Q_sorted = Q[idx]
Dcum = similar(D.density)
Dcum[idx] .= 100cumsum(Q_sorted)
ux, D, Dcum

cmap = cgrad(:oslo, categorical=true, rev=true)

contourf(D.x * ux .|> mmol/m^3, D.y * ux .|> mmol/m^3, Dcum,
         color=cmap, lc=:black, la=0, lw=0, levels=0:5:100)
plot!(collect(boundary)*ux, collect(boundary)*ux, c=:red, label="1:1")
plot!(collect(boundary)*ux, collect(boundary)*ux*0.9, c=:red, linestyle=:dot, label="±10%")
plot!(collect(boundary)*ux, collect(boundary)*ux*1.1, c=:red, linestyle=:dot, label="", 
      xlab="Observed Si(OH)₄", ylab="Modelled Si(OH)₄", legend=:bottomright)
xlims!(boundary)
ylims!(boundary)

# Add label






# Plot bSi distribution
WOA_blue   = RGB(0.231, 0.631, 0.922)
WOA_red    = RGB(0.957, 0.271, 0.188)
WOA_yellow = RGB(0.984 , 0.929, 0.298)
cmap = cgrad([WOA_blue, WOA_yellow, WOA_red])
plotverticalmean(PSi * mol/m^3 .|> μmol/m^3, grd, c=cmap, title="vertical mean opal", clim=(0,10))
=#



# Save bSi distribution



