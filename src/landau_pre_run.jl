using Statistics, Dates, JLD2, FileIO, IterTools
using Formatting, Printf

run_file = joinpath(@__DIR__, "../Run.jl")
julia_bin = "/usr/people/snirgaz/omerp/.juliaup/bin/julia"

######################### Generate slurm.txt #################################
templateSLURM = FormatExpr("""#!/bin/bash
#SBATCH --job-name={1}
#SBATCH --output={2}/%A_%a.out
#SBATCH --error={2}/%A_%a.error
#SBATCH -c 1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --requeue
#SBATCH -A snirgaz-account
#SBATCH -p snirq
#SBATCH --time=168:00:00
#SBATCH --qos=normal
#SBATCH --array=1-{3}
#SBATCH --mail-type=end

export DIR={2}
# השינוי: אנו מצביעים לקובץ sim.jld2 ספציפית
$julia_bin --threads 1 --check-bounds=no -O3 $run_file --sim_path="\$DIR/sim.jld2"
""")
###############################################################################

include("../simulation_params.jl") 

# נתיב דאטה
path = "/usr/people/snirgaz/omerp/data/sim_data/CSS_spin_rep/" * geo_name

# טעינת גיאומטריה מהתיקייה היחסית בשרת (../graphs)
geo_file_path = joinpath(@__DIR__, "../graphs", geo_name * ".jld2")
if !isfile(geo_file_path)
    error("Geometry file not found at: $geo_file_path")
end
loaded_geo = load(geo_file_path)
geo_dict = haskey(loaded_geo, "geometry_dict") ? loaded_geo["geometry_dict"] : loaded_geo

version_tracker = Dict{String, Int}()

for (β, h, J, ϵ, n_meas, n_sweep, n_therm) in IterTools.product(betas, hs, Jz, epsilons, nm_meas, nm_sweep, nm_therm)
    
    base_name_params = @sprintf("beta_%.2f_h_%.2f_eps_%.4f", β, h, ϵ)
    
    # ניהול גרסאות
    if haskey(version_tracker, base_name_params)
        version_tracker[base_name_params] += 1
        version = version_tracker[base_name_params]
    else
        version = 1
        c_path_check = joinpath(path, "$(base_name_params)_v$(version)")
        while isdir(c_path_check)
            version += 1
            c_path_check = joinpath(path, "$(base_name_params)_v$(version)")
        end
        version_tracker[base_name_params] = version
    end

    c_path = joinpath(path, "$(base_name_params)_v$(version)")
    mkpath(c_path)
    
    # שמירת נתוני ההרצה
    c_sim_data = Dict(
        "geometry" => geo_dict, 
        "beta"     => β,
        "h"        => h,
        "Jz"       => J,
        "epsilon"  => ϵ,
        "nm_meas"  => n_meas,
        "nm_sweep" => n_sweep,
        "nm_therm" => n_therm
    )
    
    sim_file = joinpath(c_path, "sim.jld2")
    save(sim_file, c_sim_data)
    
    slurm_file_name = joinpath(c_path, "slurm.txt")
    slurm_job_name = @sprintf("CSS_%s_b%.1f_v%d", geo_name, β, version)
    
    open(slurm_file_name, "w") do slurm_file
        printfmt(slurm_file, templateSLURM, slurm_job_name, c_path, 1)
    end

    println("🚀 Submitting job: $slurm_job_name")
    run(`sbatch $slurm_file_name`)
end