using JLD2
using Plots
using Statistics
using Printf

# ==========================================
# הגדרות נתיבים
# ==========================================
const SIM_DIR = "/Users/omerp/CSS_spin_rep/ED/two_cubes"
const ED_FILE = joinpath(@__DIR__, "two_cubes_ED.jld2")
const OUTPUT_DIR = joinpath(SIM_DIR, "plots")

mkpath(OUTPUT_DIR)

# ==========================================
# פונקציות עזר
# ==========================================

function parse_sim_filename(filename)
    m = match(r"res_beta_([\d\.]+)_eps_inv_(\d+)_h_([\d\.]+)\.jld2", filename)
    if m === nothing
        return nothing
    end
    beta = parse(Float64, m.captures[1])
    eps_inv = parse(Int, m.captures[2])
    h = parse(Float64, m.captures[3])
    return (beta, eps_inv, h)
end

function plot_comparison(beta, h, exact_flux, exact_corr, matching_files)
    println("  Plotting for beta=$beta, h=$h (found $(length(matching_files)) epsilon files)...")

    # 1. Flux Average (Linear X, Linear Y) - נשאר רגיל כדי לראות את הערך
    p1 = plot(title="Flux Avg (β=$beta, h=$h)", 
              ylabel="<Flux>", 
              legend=:outertopright, grid=true)
    hline!(p1, [exact_flux], label="ED Exact", color=:black, linestyle=:dash, lw=2)

    # 2. Flux Error (Log-Log) <--- שינוי כאן
    p2 = plot(title="Flux Error",
              ylabel="log(|ΔFlux|)", xlabel="log(MC Steps)",
              yaxis=:log, xaxis=:log, # גם ציר X לוגריתמי
              legend=:outertopright, grid=true)

    # 3. Correlation Average (Linear X, Linear Y)
    p3 = plot(title="Corr Avg (β=$beta, h=$h)", 
              ylabel="<Corr>", xlabel="MC Steps",
              legend=:outertopright, grid=true)
    hline!(p3, [exact_corr], label="ED Exact", color=:black, linestyle=:dash, lw=2)

    # 4. Correlation Error (Log-Log) <--- שינוי כאן
    p4 = plot(title="Corr Error",
              ylabel="log(|ΔCorr|)", xlabel="log(MC Steps)",
              yaxis=:log, xaxis=:log, # גם ציר X לוגריתמי
              legend=:outertopright, grid=true)


    # === טעינת נתונים וציור ===
    sort!(matching_files, by = x -> x[2]) 
    colors = theme_palette(:auto) 
    max_len = 0

    for (i, (file_path, eps_inv)) in enumerate(matching_files)
        data = load(file_path)
        meas = data["measurements"]
        
        # Flux
        flux_samples = meas.avg_spatial_flux
        steps = 1:length(flux_samples)
        max_len = max(max_len, length(steps))
        
        run_avg_flux = cumsum(flux_samples) ./ steps
        res_flux = abs.(run_avg_flux .- exact_flux) .+ 1e-20 

        # Correlation
        corr_samples = meas.flux_corr
        steps_c = 1:length(corr_samples)
        run_avg_corr = cumsum(corr_samples) ./ steps_c
        res_corr = abs.(run_avg_corr .- exact_corr) .+ 1e-20

        lbl = "ε = 1/$eps_inv"
        c = colors[mod1(i, length(colors))]

        plot!(p1, steps, run_avg_flux, label=lbl, lw=1.5, color=c, alpha=0.8)
        plot!(p2, steps, res_flux,     label=lbl, lw=1,   color=c, alpha=0.6)
        
        plot!(p3, steps_c, run_avg_corr, label=lbl, lw=1.5, color=c, alpha=0.8)
        plot!(p4, steps_c, res_corr,     label=lbl, lw=1,   color=c, alpha=0.6)
    end

    # Theoretical line 1/sqrt(N)
    if max_len > 0
        xs = 100:max_len 
        theory = 0.5 ./ sqrt.(xs) 
        
        # הקו התיאורטי ייראה ישר בגרף log-log
        plot!(p2, xs, theory, label="~ 1/√N", color=:black, linestyle=:dot, lw=2)
        plot!(p4, xs, theory, label="~ 1/√N", color=:black, linestyle=:dot, lw=2)
    end

    l = @layout [a b; c d]
    p_final = plot(p1, p2, p3, p4, layout=l, size=(1200, 900), margin=8Plots.mm)
    
    out_name = @sprintf("conv_beta_%.1f_h_%.2f.png", beta, h)
    save_path = joinpath(OUTPUT_DIR, out_name)
    savefig(p_final, save_path)
end

function main()
    println("Scanning simulation directory: $SIM_DIR")
    if !isdir(SIM_DIR)
        error("Simulation dir not found!")
    end

    # מיפוי קבצים
    sim_registry = []
    for f in readdir(SIM_DIR)
        if endswith(f, ".jld2")
            parsed = parse_sim_filename(f)
            if parsed !== nothing
                push!(sim_registry, (joinpath(SIM_DIR, f), parsed...))
            end
        end
    end
    println("Found $(length(sim_registry)) valid simulation files.")

    println("Loading ED data from: $ED_FILE")
    if !isfile(ED_FILE)
        error("ED file not found!")
    end
    ed_data = load(ED_FILE)
    
    hs_vector = ed_data["hs"]

    # לולאה על ה-ED keys
    for k in keys(ed_data)
        if k == "hs" continue end
        
        beta_val = tryparse(Float64, k)
        if beta_val === nothing continue end

        beta_group = ed_data[k]
        
        # בדיקה שיש את שני הנתונים
        if !haskey(beta_group, "F_avg") || !haskey(beta_group, "F_corr")
            println("Missing F_avg or F_corr for beta=$beta_val")
            continue
        end
        
        f_avgs = beta_group["F_avg"]
        f_corrs = beta_group["F_corr"] # <--- שליפת הקורלציות מה-ED

        for (i, h_val) in enumerate(hs_vector)
            exact_flux = f_avgs[i]
            exact_corr = f_corrs[i] # <--- הערך המדויק לקורלציה
            
            # חיפוש קבצים מתאימים
            relevant_files = []
            for (path, s_beta, s_eps_inv, s_h) in sim_registry
                if isapprox(s_beta, beta_val, atol=1e-5) && isapprox(s_h, h_val, atol=1e-5)
                    push!(relevant_files, (path, s_eps_inv))
                end
            end

            if isempty(relevant_files)
                continue
            end

            plot_comparison(beta_val, h_val, exact_flux, exact_corr, relevant_files)
        end
    end
    
    println("\nDone! Plots saved to $OUTPUT_DIR")
end

main()