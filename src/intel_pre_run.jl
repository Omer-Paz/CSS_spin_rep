using Statistics, Dates, JLD2, FileIO, IterTools
using Formatting, Printf

# --- שינוי 1: תיקון הנתיב ל-Run.jl (כמו ב-Landau) ---
run_file = joinpath(@__DIR__, "Run.jl")

# נתיב הג'וליה הספציפי של Intel (נשמר)
julia_bin = "/homes/omerp/.juliaup/bin/julia"

######################### Generate slurm.txt #################################
templateSLURM = FormatExpr("""#!/bin/bash
#SBATCH --job-name={1}
#SBATCH --output={2}/%A_%a.out
#SBATCH --error={2}/%A_%a.error
#SBATCH -c 1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --requeue
#SBATCH --time=168:00:00
#SBATCH --qos=normal
#SBATCH --array=1-{3}
#SBATCH --mail-type=end
#SBATCH --mail-user=omer.paz2@mail.huji.ac.il

export DIR={2}
$julia_bin --threads 1 --check-bounds=no -O3 $run_file --sim_path="\$DIR/sim.jld2"
""")
###############################################################################

# --- שינוי 2: טעינת הפרמטרים מהתיקייה הנוכחית (src) ---
include("simulation_params.jl") 

# נתיב שמירת הנתונים ב-Intel (נשמר המיקום המקורי Z2_Gauge)
path = joinpath("/homes/omerp/sim_data/CSS_spin_rep/", geo_folder, geo_name)

# --- שינוי 3: טעינת גיאומטריה מתיקיית graphs (במקום geometries) ---
geo_file_path = joinpath(dirname(@__DIR__),geo_folder, geo_name * ".jld2")
@show geo_file_path
if !isfile(geo_file_path)
    error("Geometry file not found at: $geo_file_path")
end
loaded_geo = load(geo_file_path)
# הוספת הבדיקה אם המילון מקונן (כמו ב-Landau)
geo_dict = haskey(loaded_geo, "geometry_dict") ? loaded_geo["geometry_dict"] : loaded_geo

# --- מעקב גרסאות ---
version_tracker = Dict{String, Int}()

# לולאה על כל הפרמטרים
# הערה: n_sweep כאן הוא הפקטור (למשל 10) שמגיע מ-simulation_params
for (β, h, J, ϵ, n_meas, n_sweep_factor, n_therm) in IterTools.product(betas, hs, Jz, epsilons, nm_meas, nm_sweep_factor, nm_therm)
    
    base_name_params = @sprintf("beta_%.2f_h_%.2f_eps_%.4f", β, h, ϵ)
    
    # 1. לוגיקה למציאת גרסה
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

    # יצירת נתיב
    c_path = joinpath(path, "$(base_name_params)_v$(version)")
    mkpath(c_path)
    
    # --- שינוי 4: חישוב מספר הצעדים הכולל לפי נפח המערכת ---
    system_vol = geo_dict["N_vertices"] + geo_dict["N_edges"]
    
    c_sim_data = Dict(
        "geometry" => geo_dict, 
        "beta"     => β,
        "h"        => h,
        "Jz"       => J,
        "epsilon"  => ϵ,
        "nm_meas"  => n_meas,
        "nm_sweep" => system_vol * n_sweep_factor, # הכפלת הנפח בפקטור
        "nm_therm" => n_therm
    )
    
    # שמירת קובץ הקלט לסימולציה
    sim_file = joinpath(c_path, "sim.jld2")
    save(sim_file, c_sim_data)
    
    # יצירת קובץ Slurm
    slurm_file_name = joinpath(c_path, "slurm.txt")
    slurm_job_name = @sprintf("CSS_%s_%s_b%.1f_v%d", geo_folder,geo_name, β, version)
    
    open(slurm_file_name, "w") do slurm_file
        printfmt(slurm_file, templateSLURM, slurm_job_name, c_path, 1)
    end

    println("🚀 Submitting job: $slurm_job_name (Folder: $geo_folder)")
    run(`sbatch $slurm_file_name`)
end