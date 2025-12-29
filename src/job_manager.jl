using Dates

# ==============================================================================
# 1. הגדרות שרת (Server Selection)
# ==============================================================================
target_server = "intel"  # אפשרויות: "intel" או "landau"

# ==============================================================================
# 2. הגדרת פרמטרים לסימולציה
# ==============================================================================
# שם קובץ הגיאומטריה (ללא סיומת) שנמצא בתיקיית geometries
geo_name = "psl_2_4" 

betas = [2.0]#[2.0,4.0,8.0,16.0] 
hs = [0.5]#collect(0.1:0.1:1.0)
Jz = [1.0]

epsilons = [1/16]#[1/8,1/16,1/32,1/64] # דיסקרטיזציה של הזמן

# הגדרות ריצה
nm_meas = [1]#[2^13]    # מספר מדידות
nm_sweep_factor = [1]     # כמה צעדים בין מדידה למדידה
nm_therm = [1]   # צעדי תרמליזציה

# תוכן הקובץ שייכתב וישלח לשרת
params_content = """
# Auto-generated parameters file for Z2 Gauge Theory
geo_name = "$geo_name"
betas = $betas
hs = $hs
Jz = $Jz
epsilons = $epsilons
nm_meas = $nm_meas
nm_sweep_factor = $nm_sweep_factor
nm_therm = $nm_therm
"""

# ==============================================================================
# 3. עדכון קבצים וגיט (Local Git Operations)
# ==============================================================================
params_path = joinpath(@__DIR__, "simulation_params.jl")
current_script_path = @__FILE__ 

write(params_path, params_content)
println("✅ Updated parameters file: $params_path")
# ...
println("🚀 Pushing code to Git...")
try
    # הוספת קובץ הפרמטרים
    run(`git add $params_path`)
    
    # הוספת כל הקבצים בתיקיית src (כולל landau_pre_run.jl המתוקן)
    src_dir = joinpath(@__DIR__, "src")
    run(`git add $src_dir`)
    
    run(`git commit -m "Auto-update params and source code"`)
    run(`git push origin main`) 
catch e
    println("⚠️ Git warning: ", e)
end
# ...# ==============================================================================
# 4. הכנת הפקודה לשרת (Dynamic Server Configuration)
# ==============================================================================
# עדכן כאן את שם התיקייה בשרת שבה יושב הפרויקט הזה
project_folder_name = "CSS_spin_rep"

if target_server == "intel"
    remote_host = "intel"
    remote_dir = "~/scripts/$project_folder_name" 
    script_to_run = "src/intel_pre_run.jl"
    
elseif target_server == "landau"
    remote_host = "landau"
    remote_dir = "~/data/scripts/$project_folder_name"
    script_to_run = "src/landau_pre_run.jl"
    
else
    error("❌ Unknown server selected: $target_server")
end

println("🔗 Connecting via Office -> $remote_host...")
println("📂 Working directory: $remote_dir")
println("📜 Script to run: $script_to_run")

remote_julia = "~/.juliaup/bin/julia"

# הפקודה שתרוץ בתוך השרת הסופי
cmd_run = "cd $remote_dir && git pull origin main && $remote_julia $script_to_run"

# הפקודה שתרוץ ב-Office (מקפצה)
cmd_office = "ssh -A -t omerp@$remote_host \"$cmd_run\""

# ==============================================================================
# 5. ביצוע ההרצה
# ==============================================================================
try
    run(`ssh -A -t office $cmd_office`)
    println("🎉 Job submitted successfully to $target_server!")
catch e
    println("❌ SSH Error: ", e)
end