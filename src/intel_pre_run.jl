using Statistics, Dates, JLD2, FileIO, IterTools
using Formatting, Printf

run_file = joinpath(@__DIR__, "../Run.jl") # בהנחה ש-Run נמצא רמה אחת מעל src או באותה תיקייה
# אם Run.jl נמצא יחד עם הקובץ הזה באותה תיקייה: joinpath(@__DIR__, "Run.jl")

julia_bin = "/homes/omerp/.juliaup/bin/julia"

######################### Generate slurm.txt #################################
templateSLURM = FormatExpr("""#!/bin/bash
#SBATCH --job-name={1}
#SBATCH --output={2}/%A_%a.out
#SBATCH --error={2}/%A_%a.error
#SBATCH -c 1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --requeue
#SBATCH --time=168:00:00
#SBATCH --qos=normal
#SBATCH --array=1-{3}
#SBATCH --mail-type=end
#SBATCH --mail-user=omer.paz2@mail.huji.ac.il

export DIR={2}
$julia_bin --threads 1 --check-bounds=no -O3 $run_file --sim_path=\$DIR
""")
###############################################################################

# טעינת הפרמטרים שנוצרו ע"י job_manager
include("../simulation_params.jl") 

# נתיב בסיס לשמירת התוצאות
path = "/homes/omerp/sim_data/Z2_Gauge/" * geo_name

# טעינת הגיאומטריה (Load Geometry Once)
# הנחה: קבצי הגיאומטריה נמצאים בתיקיית geometries יחסית לסקריפט או בנתיב קבוע
geo_file_path = joinpath(@__DIR__, "../geometries", geo_name * ".jld2")
if !isfile(geo_file_path)
    error("Geometry file not found at: $geo_file_path")
end
loaded_geo = load(geo_file_path)
# אם הקובץ מכיל מפתח ראשי, שלוף אותו. אם הוא המילון עצמו:
geo_dict = haskey(loaded_geo, "geometry_dict") ? loaded_geo["geometry_dict"] : loaded_geo

# --- מעקב גרסאות ---
version_tracker = Dict{String, Int}()

# לולאה על כל הפרמטרים
for (β, h, J, ϵ, n_meas, n_sweep, n_therm) in IterTools.product(betas, hs, Jz, epsilons, nm_meas, nm_sweep, nm_therm)
    
    # שם התיקייה הבסיסי
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
    
    # הכנת המילון עבור Run.jl (חייב להתאים למפתחות ב-Run.jl)
    c_sim_data = Dict(
        "geometry" => geo_dict, # שים לב: זה המפתח ש-Run.jl מחפש
        "beta"     => β,
        "h"        => h,
        "Jz"       => J,
        "epsilon"  => ϵ,
        "nm_meas"  => n_meas,
        "nm_sweep" => n_sweep,
        "nm_therm" => n_therm
    )
    
    # שמירת קובץ הקלט לסימולציה
    sim_file = joinpath(c_path, "sim.jld2")
    save(sim_file, c_sim_data)
    
    # יצירת קובץ Slurm
    slurm_file_name = joinpath(c_path, "slurm.txt")
    slurm_job_name = @sprintf("Z2_%s_b%.1f_h%.2f_v%d", geo_name, β, h, version)
    
    open(slurm_file_name, "w") do slurm_file
        printfmt(slurm_file, templateSLURM, slurm_job_name, c_path, 1)
    end

    println("🚀 Submitting job: $slurm_job_name")
    run(`sbatch $slurm_file_name`)
end