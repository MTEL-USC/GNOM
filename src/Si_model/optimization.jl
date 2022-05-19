# solve it once
sol = solve(prob, CTKAlg(), preprint="Nd & εNd solve ", τstop=ustrip(s, 1e3Myr))

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
τstop = ustrip(s, 1e3Myr)
mem = F1Method.initialize_mem(F, ∇ₓf, ∇ₓF, x0, λ, CTKAlg(); preprint="mem ", τstop=τstop)

function objective(λ)
    p = λ2p(SiParams, λ) ; @show p
    F1Method.objective(f, F, ∇ₓF, mem, λ, CTKAlg(), preprint="obj ", τstop=τstop)
end
gradient(λ) = F1Method.gradient(f, F, ∇ₓf, ∇ₓF, mem, λ, CTKAlg(), preprint="grad", τstop=τstop)
hessian(λ) = F1Method.hessian(f, F, ∇ₓf, ∇ₓF, mem, λ, CTKAlg(), preprint="hess ", τstop=τstop)

# Reduced g_tol for optimization
opt = Optim.Options(store_trace=false, show_trace=true, extended_trace=false, g_tol=1e-5)

#================================================
Optimization run
================================================#
results = optimize(objective, gradient, hessian, λ, NewtonTrustRegion(), opt; inplace=false)

p_optimized = λ2p(SiParams, results.minimizer)
prob_optimized = SteadyStateProblem(F, x0, p_optimized)
s_optimized = solve(prob_optimized, CTKAlg(), τstop=ustrip(s, 1e3Myr)).u

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



