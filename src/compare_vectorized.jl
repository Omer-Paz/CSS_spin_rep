using JLD2
using Statistics
using Printf
using Random

# טעינת הקבצים שלך
include("Config.jl") 
include("Measurements.jl") 
include("Params.jl") 
include("Metropolis.jl") 

# --- פונקציות עזר (זהות למה שהיה) ---

function load_geometry(path::String)
    if !isfile(path)
        error("Geometry file not found at: $path")
    end
    println("Loading geometry from $path...")
    return load(path)
end

function initialize_simulation(beta, h, epsilon, Jz, nm_therm, nm_meas, nm_sweep, geo_dict)
    M = Int(round(beta / epsilon))
    params = SimParams(beta, M, Jz, h, epsilon, nm_sweep, nm_therm, nm_meas)
    conf = SimConfig(geo_dict, params.M)
    N_vert = length(conf.vert_to_edge)
    meas = init_measurements(params.nm_meas, N_vert)
    
    # הדפסה קצרה
    @printf("Init: beta=%.1f, eps=1/%d (M=%d), h=%.2f\n", params.beta, Int(1/epsilon), params.M, params.h)
    return params, conf, meas
end

function run_loop!(meas::Measurements, conf::SimConfig, params::SimParams, move_func!::Function)
    # פונקציה כללית שמריצה תרמליזציה + מדידות לצורך הבדיקה הפיזיקלית
    # 1. Thermalization
    for _ in 1:params.nm_therm
        for _ in 1:params.nm_sweep
            move_func!(conf, params)
        end
    end
    # 2. Measurements
    for i in 1:params.nm_meas
        for _ in 1:params.nm_sweep
            move_func!(conf, params)
        end
        measure_all!(meas, conf)
    end
end

function analyze_results(meas::Measurements)
    avg_flux = mean(meas.avg_spatial_flux)
    err_flux = std(meas.avg_spatial_flux) / sqrt(length(meas.avg_spatial_flux))
    avg_corr = mean(meas.flux_corr)
    err_corr = std(meas.flux_corr) / sqrt(length(meas.flux_corr))
    return avg_flux, err_flux, avg_corr, err_corr
end

# --- פונקציות ההשוואה ---

function compare_moves_performance(conf::SimConfig, params::SimParams, n_sweeps_bench::Int=500)
    println("\n=== ⏱️  Performance Benchmark ($n_sweeps_bench sweeps) ===")
    
    # 1. Warmup (Compilation)
    print("  Warming up... ")
    sweep_move!(conf, params)
    sweep_move_vectorized!(conf, params)
    println("Done.")

    # 2. Sequential Benchmark
    println("  Running Sequential move...")
    t_seq = @elapsed for _ in 1:n_sweeps_bench
        sweep_move!(conf, params)
    end
    speed_seq = n_sweeps_bench / t_seq

    # 3. Vectorized Benchmark
    println("  Running Vectorized move...")
    t_vec = @elapsed for _ in 1:n_sweeps_bench
        sweep_move_vectorized!(conf, params)
    end
    speed_vec = n_sweeps_bench / t_vec

    println("  ------------------------------------------------")
    @printf("  Sequential Time: %.4f s (%.1f sweeps/sec)\n", t_seq, speed_seq)
    @printf("  Vectorized Time: %.4f s (%.1f sweeps/sec)\n", t_vec, speed_vec)
    println("  ------------------------------------------------")
    @printf("  🚀 Speedup Factor: x%.2f\n", speed_seq / speed_vec)
    println("  ------------------------------------------------")
end

function compare_moves_physics(beta, h, epsilon, Jz, geo_dict)
    println("\n=== 🔬 Physics Validation Check ===")
    
    # פרמטרים לריצה קצרה אך מספקת לסטטיסטיקה בסיסית
    short_therm = 500
    short_meas = 2000
    short_sweep = 10
    
    println("  Params: Therm=$short_therm, Meas=$short_meas, Sweep=$short_sweep")

    # --- הרצה 1: רגיל (Sequential) ---
    println("  1. Running Sequential Simulation...")
    p1, c1, m1 = initialize_simulation(beta, h, epsilon, Jz, short_therm, short_meas, short_sweep, geo_dict)
    @time run_loop!(m1, c1, p1, sweep_move!)
    f1, ef1, c_val1, ec1 = analyze_results(m1)
    
    # --- הרצה 2: וקטורי (Vectorized) ---
    println("  2. Running Vectorized Simulation...")
    p2, c2, m2 = initialize_simulation(beta, h, epsilon, Jz, short_therm, short_meas, short_sweep, geo_dict)
    @time run_loop!(m2, c2, p2, sweep_move_vectorized!)
    f2, ef2, c_val2, ec2 = analyze_results(m2)

    # --- הדפסת טבלה ---
    println("\n  -------------------------------------------------------------")
    println("  Variable      | Sequential            | Vectorized")
    println("  -------------------------------------------------------------")
    @printf("  Flux          | %.5f ± %.5f   | %.5f ± %.5f\n", f1, ef1, f2, ef2)
    @printf("  Correlation   | %.5f ± %.5f   | %.5f ± %.5f\n", c_val1, ec1, c_val2, ec2)
    println("  -------------------------------------------------------------")
    
    # בדיקת תאימות (מרחק של 3 סטיות תקן משוקללות)
    sigma_diff = sqrt(ef1^2 + ef2^2)
    diff = abs(f1 - f2)
    
    if diff < 3 * sigma_diff
        println("  ✅ SUCCESS: Physics results match (difference < 3σ).")
    else
        println("  ⚠️ WARNING: Physics mismatch! Difference is $(round(diff/sigma_diff, digits=1))σ.")
    end
end

function main()
    # === הגדרת פרמטרים לבדיקה ===
    # בחר פרמטרים "כבדים" מספיק כדי להרגיש את ההבדל, אבל לא כבדים מדי
    beta = 2.0
    epsilon = 1/16
    h = 0.5
    Jz = 1.0
    
    # נתיבים
    geo_path = joinpath(dirname(@__DIR__),"graphs", "psl_2_4.jld2")
    geo_dict = load_geometry(geo_path)

    # 1. יצירת סביבה לבדיקת מהירות
    # אנחנו מאתחלים עם מעט צעדים רק בשביל ליצור את ה-Config
    println("Initializing configuration for benchmarks...")
    params, conf, _ = initialize_simulation(beta, h, epsilon, Jz, 10, 10, 10, geo_dict)
    
    # 2. הרצת בדיקת המהירות
    compare_moves_performance(conf, params, 1000) # מספר ה-sweeps למדידה
    
    # 3. הרצת בדיקת הפיזיקה
    compare_moves_physics(beta, h, epsilon, Jz, geo_dict)
    
    println("\nDone.")
end

main()