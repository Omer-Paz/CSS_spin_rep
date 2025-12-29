using JLD2
using Statistics
using Printf

include("Config.jl") 
include("Measurements.jl") 
include("Params.jl") 
include("Metropolis.jl") 

function load_geometry(path::String)
    if !isfile(path)
        error("Geometry file not found at: $path")
    end
    println("Loading geometry from $path...")
    return load(path)
end

# שינוי קטן: הפונקציה מקבלת כעת את geo_dict הטעון במקום נתיב
function initialize_simulation(beta, h, epsilon, Jz, nm_therm, nm_meas, nm_sweep, geo_dict)
    M = Int(round(beta / epsilon))
    params = SimParams(beta, M, Jz, h, epsilon, nm_sweep, nm_therm, nm_meas)
    conf = SimConfig(geo_dict, params.M)
    meas = init_measurements(params.nm_meas)
    
    # הדפסה קצרה יותר כדי לא להציף את הלוג בריצה ארוכה
    @printf("Init: beta=%.1f, eps=1/%d (M=%d), h=%.2f\n", params.beta, Int(1/epsilon), params.M, params.h)
    return params, conf, meas
end

function run_thermalization!(conf::SimConfig, params::SimParams)
    println("  Thermalizing...") 
    for _ in 1:params.nm_therm
        sweep_move!(conf, params)
    end
end

function run_measurement_loop!(meas::Measurements, conf::SimConfig, params::SimParams)
    println("  Measuring...")
    for i in 1:params.nm_meas
        for _ in 1:params.nm_sweep
            sweep_move!(conf, params)
        end
        measure_all!(meas, conf)
    end
end

function analyze_results(meas::Measurements)
    avg_flux = mean(meas.avg_spatial_flux)
    err_flux = std(meas.avg_spatial_flux) / sqrt(length(meas.avg_spatial_flux))
    avg_corr = mean(meas.flux_corr)
    err_corr = std(meas.flux_corr) / sqrt(length(meas.flux_corr))
    
    @printf("  >> Result: Flux=%.4f±%.4f, Corr=%.4f±%.4f\n", avg_flux, err_flux, avg_corr, err_corr)
end

function main()
    # === הגדרת פרמטרים לסריקה ===
    betas = [2.0]
    epsilons = [1/16]
    hs = [0.5]#collect(0.1:0.1:1.0)
    
    Jz = 1.0
    nm_therm = 2^10
    nm_meas = 2^15
    nm_sweep = 400
    
    # === הגדרת נתיבים ===
    geo_path = joinpath(dirname(@__DIR__),"graphs", "psl_2_4.jld2")
    output_dir = dirname(@__DIR__)#"/Users/omerp/CSS_spin_rep/ED/two_cubes"
    
    # יצירת התיקייה אם אינה קיימת
    mkpath(output_dir)
    
    # טעינת הגאומטריה פעם אחת בלבד
    geo_dict = load_geometry(geo_path)
    
    total_runs = length(betas) * length(epsilons) * length(hs)
    curr_run = 0

    println("Starting Sweep. Total runs: $total_runs")
    println("Output directory: $output_dir")

    for beta in betas
        for epsilon in epsilons
            for h in hs
                curr_run += 1
                println("\nRun $curr_run / $total_runs:")
                params, conf, meas = initialize_simulation(beta, h, epsilon, Jz, nm_therm, nm_meas, nm_sweep, geo_dict)
                run_thermalization!(conf, params)    
                run_measurement_loop!(meas, conf, params)    
                analyze_results(meas)
                filename = @sprintf("res_beta_%.1f_eps_inv_%d_h_%.2f_4.jld2", beta, Int(round(1/epsilon)), h)
                save_path = joinpath(output_dir, filename)
                save(save_path, Dict("measurements" => meas, "params" => params))
            end
        end
    end
    println("\nAll runs completed.")
end
main()