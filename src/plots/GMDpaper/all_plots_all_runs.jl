use_GLMakie = false
run_num = 1
resetup = false
reload = false
retools = false
include("../plots_setup_Nd.jl")

run_nums = [parse(Int, match(r"run(\d+)_(.*).jld2", f).captures[1]) for f in readdir(archive_path) if contains(f, "jld2")]
for i in run_nums
    println("==========================================")
    println("==========================================\n")
    global run_num = i # necessary for the include statement to see run_num?
    println("Loading run $run_num...")
    reload = true # force a reload of the specific run $num_run
    include("../plots_setup_Nd.jl")
    include("GMD_diagnostics_setup.jl") # force re computations of diagnostics
    println("Run $run_num Params:")
    @show p
    include("all_plots.jl")
    println("\n==========================================")
    println("==========================================")
end


