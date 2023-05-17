# solve it once
sol = solve(prob, CTKAlg(), preprint="Nd & εNd solve ", τstop=ustrip(s, 1e3Myr))

#================================================
Optimization setup
================================================#

include("obs.jl")

# isotope function returns [Nd] and εNd from χNd and χRNd
modify(DNd, ⁱDNd) = (DNd, ⁱDNd ./ DNd ./ R_CHUR .- 1)

# Weights for objective function
const ωs = (1.0, 1.0) # the weight for the mismatches (note the new API reduced to observed parts)
const ωp = 1e-4
const obs = (DNdobs, εNdobs)

f, ∇ₓf = f_and_∇ₓf(ωs, ωp, grd, modify, obs, Params)

# Randomized starting parameter from distributions
p = Params(; ((k=>rand(d)) for (k,d) in zip(fieldnames(Params),prior(Params)) if !isnothing(d))...)

# create directories to write in
mkpath(output_path)
mkpath(data_path)
mkpath(archive_path)
# Check previous runs and get new run number
run_num = let
    previous_run_nums = [parse(Int, match(r"run(\d+)", f).captures[1]) for f in readdir(archive_path) if contains(f, "run")]
    run_num = 1
    while run_num ∈ previous_run_nums
        run_num += 1
    end
    run_num
end
println("==========================================")
println("This is run $run_num with starting Params:")
touch(joinpath(archive_path, "run$run_num")) # touch an empty file to make sure run number is unique
@show p
println("==========================================")

# Use F1 for gradient and Hessian
λ = p2λ(p)
τstop = ustrip(s, 1e3Myr)
mem = F1Method.initialize_mem(F, ∇ₓf, x, λ, CTKAlg(); preprint="mem ", τstop=τstop)

function objective(λ)
    p = λ2p(Params, λ) ; @show p
    F1Method.objective(f, F, mem, λ, CTKAlg(), preprint="obj ", τstop=τstop)
end
gradient(λ) = F1Method.gradient(f, F, ∇ₓf, mem, λ, CTKAlg(), preprint="grad", τstop=τstop)
hessian(λ) = F1Method.hessian(f, F, ∇ₓf, mem, λ, CTKAlg(), preprint="hess ", τstop=τstop)

# Reduced g_tol for optimization
opt = Optim.Options(store_trace=false, show_trace=true, extended_trace=false, g_tol=1e-4)

#================================================
Optimization run
================================================#
results = optimize(objective, gradient, hessian, λ, NewtonTrustRegion(), opt; inplace=false)

p_optimized = λ2p(Params, results.minimizer)
prob_optimized = SteadyStateProblem(F, x, p_optimized)
s_optimized = solve(prob_optimized, CTKAlg(), τstop=ustrip(s, 1e3Myr)).u

# TODO find a more generic approach to save this data... Maybe I can use datadeps?

DNd, ⁱDNd = unpack_tracers(s_optimized, grd)
DNd, εNd = modify(DNd, ⁱDNd)
tp_opt = AIBECS.table(p_optimized)



jldsave(joinpath(archive_path, "optimized_output_$(headcommit)_run$(run_num)_$(circname).jld2"); DNd, εNd, εNdobs, DNdobs, s_optimized, tp_opt)

println("Done!")

# print optimized parameters in md table
open(joinpath(archive_path, "optimized_parameters_$(headcommit)_run$(run_num)_$(circname).md"), "w") do f
    pretty_table(f, AIBECS.table(p_optimized)[:,[[1,2,3,4];6:end]], nosubheader=true, tf=tf_markdown)
end


