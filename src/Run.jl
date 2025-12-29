using JLD2
using Statistics
using Printf
using ArgParse
include("Config.jl") 
include("Measurements.jl") 
include("Params.jl") 
include("Metropolis.jl") 


function get_measurement_filename(output_dir::String)
    base_name = "measurements.jld2"
    full_path = joinpath(output_dir, base_name)
    if !isfile(full_path)
        return base_name
    end
    i = 1
    while isfile(joinpath(output_dir, "measurements_$(i).jld2"))
        i += 1
    end
    return "measurements_$(i).jld2"
end

function save_checkpoint(output_dir::String, meas::Measurements, conf::SimConfig, params::SimParams, meas_filename::String)
    jldsave(joinpath(output_dir, meas_filename); measurements=meas)
    jldsave(joinpath(output_dir, "config.jld2"); config=conf)
    params_path = joinpath(output_dir, "params.jld2")
    if !isfile(params_path)
        jldsave(params_path; params=params)
    end
end

function run_therm!(conf::SimConfig, params::SimParams)
    for i in 1:params.nm_therm
        for _=1:params.nm_sweep
            sweep_move_vectorized!(conf, params)
        end
    end
end

function run_measurements!(meas::Measurements, conf::SimConfig, params::SimParams,output_dir::String, meas_filename::String)
    for i in 1:params.nm_meas
        for _ in 1:params.nm_sweep
            sweep_move_vectorized!(conf, params)
        end
        measure_all!(meas, conf)
        if i % 1000 == 0
            @printf("[Checkpoint] %d/%d measurements saved.\n", i, params.nm_meas)
            flush(stdout)
            save_checkpoint(output_dir, meas, conf, params, meas_filename)
        end
    end
    save_checkpoint(output_dir, meas, conf, params, meas_filename)
    println("✅ Finished measurements. Saved to $meas_filename")
end

function initialize_simulation(beta, h, epsilon, Jz, nm_therm, nm_meas, nm_sweep, geo_dict)
    M = Int(round(beta / epsilon))
    params = SimParams(beta, M, Jz, h, epsilon, nm_sweep, nm_therm, nm_meas)
    conf = SimConfig(geo_dict, params.M)
    N_vert = length(conf.vert_to_edge)
    meas = init_measurements(params.nm_meas,N_vert)    
    @printf("Init: beta=%.1f, eps=1/%d (M=%d), h=%.2f\n", params.beta, Int(1/epsilon), params.M, params.h)
    flush(stdout)
    return params, conf, meas
end

function initialize_new_run(sim_data::Dict)

    geo_dict = sim_data["geometry"] 
    beta     = sim_data["beta"]
    h        = sim_data["h"]
    epsilon  = sim_data["epsilon"]
    Jz       = sim_data["Jz"]
    nm_therm = sim_data["nm_therm"]
    nm_meas  = sim_data["nm_meas"]
    nm_sweep = sim_data["nm_sweep"]
    
    M = Int(round(beta / epsilon))
    
    params = SimParams(beta, M, Jz, h, epsilon, nm_sweep, nm_therm, nm_meas)
    conf = SimConfig(geo_dict, params.M)
    N_vert = length(conf.vert_to_edge)
    meas = init_measurements(params.nm_meas,N_vert)
    
    @printf("✨ Init NEW run: beta=%.2f, M=%d, h=%.2f\n", beta, M, h)
    return params, conf, meas
end

function run_sim(sim_input_path::String)
    output_dir = dirname(sim_input_path)
    config_path = joinpath(output_dir, "config.jld2")
    params_path = joinpath(output_dir, "params.jld2")
    is_continuation = isfile(config_path) && isfile(params_path)
    local params, conf, meas
    if is_continuation
        println("🔄 Found existing configuration. Resuming run...")
        loaded_conf = load(config_path, "config")
        loaded_params = load(params_path, "params")
        N_vert = length(loaded_conf.vert_to_edge)
        meas = init_measurements(loaded_params.nm_meas,N_vert)
        params = loaded_params
        conf = loaded_conf        
    else
        input_data = load(sim_input_path)
        params, conf, meas = initialize_new_run(input_data)
        jldsave(joinpath(output_dir, "params.jld2"); params=params)
        run_therm!(conf, params)
    end
    
    # קביעת שם קובץ הפלט (measurements.jld2 או measurements_1.jld2 וכו')
    current_meas_filename = get_measurement_filename(output_dir)
    
    # הרצת המדידות
    run_measurements!(meas, conf, params, output_dir, current_meas_filename)
end

s = ArgParseSettings()
@add_arg_table s begin
     "--sim_path"
         help = "Path to the sim.jld2 input file"
         required = true
end
parsed_args = parse_args(ARGS, s)
run_sim(parsed_args["sim_path"])