using JLD2
using Plots
using Statistics
using Printf

# ==========================================
# 1. הגדרות נתיבים
# ==========================================
const SIM_DIR = "/Users/omerp/CSS_spin_rep/ED/two_cubes"
const ED_FILE = joinpath(@__DIR__, "two_cubes_ED.jld2")
const OUTPUT_DIR = joinpath(SIM_DIR, "plots")

mkpath(OUTPUT_DIR)

# ==========================================
# 2. פונקציות עזר ופארסינג
# ==========================================

function parse_sim_filename(filename)
    # Regex לזיהוי הפרמטרים והסיומת
    m = match(r"res_beta_([\d\.]+)_eps_inv_(\d+)_h_([\d\.]+)(.*)\.jld2", filename)
    if m === nothing
        return nothing
    end
    
    beta = parse(Float64, m.captures[1])
    eps_inv = parse(Int, m.captures[2])
    h = parse(Float64, m.captures[3])
    remainder = m.captures[4]
    
    # זיהוי האם הקובץ הוא Vectorized
    is_vec = occursin("vectorized", remainder)
    
    return (beta, eps_inv, h, is_vec)
end

# ==========================================
# 3. פונקציות ציור
# ==========================================

# פונקציה גנרית שמציירת קבוצת קבצים (או רגיל או וקטורי) מול ה-ED
# כולל גרף שגיאה ב-LogLog
function plot_group_vs_ed(beta, h, exact_flux, exact_corr, matching_files, mode_prefix)
    if isempty(matching_files) return end
    
    # הגדרת כותרות וצבעים לפי המוד
    if mode_prefix == "vec"
        type_title = "Vectorized"
        out_prefix = "conv_vec"
    else
        type_title = "Sequential"
        out_prefix = "conv_seq"
    end
    
    println("  Creating $type_title plots for beta=$beta, h=$h...")

    # --- פאנל 1: ממוצע שטף (ליניארי) ---
    p1 = plot(title="$type_title Flux (β=$beta)", ylabel="<Flux>", 
              legend=:outertopright, grid=true, guidefontsize=8)
    hline!(p1, [exact_flux], label="ED Exact", color=:black, linestyle=:dash, lw=2)

    # --- פאנל 2: שגיאת שטף (Log-Log) ---
    p2 = plot(title="Flux Error (Log-Log)", ylabel="log(|ΔFlux|)", xlabel="log(Steps)",
              yaxis=:log, xaxis=:log, legend=:outertopright, grid=true, guidefontsize=8)

    # --- פאנל 3: ממוצע קורלציה (ליניארי) ---
    p3 = plot(title="$type_title Corr (β=$beta)", ylabel="<Corr>", xlabel="Steps",
              legend=:outertopright, grid=true, guidefontsize=8)
    hline!(p3, [exact_corr], label="ED Exact", color=:black, linestyle=:dash, lw=2)

    # --- פאנל 4: שגיאת קורלציה (Log-Log) ---
    p4 = plot(title="Corr Error (Log-Log)", ylabel="log(|ΔCorr|)", xlabel="log(Steps)",
              yaxis=:log, xaxis=:log, legend=:outertopright, grid=true, guidefontsize=8)

    # מיון הקבצים לפי גודל הצעד (אפסילון)
    sort!(matching_files, by = x -> x[2]) 
    colors = theme_palette(:auto) 
    max_len = 0

    for (i, (file_path, eps_inv)) in enumerate(matching_files)
        data = load(file_path)
        meas = data["measurements"]
        
        # --- עיבוד Flux ---
        flux_samples = meas.avg_spatial_flux
        steps = 1:length(flux_samples)
        max_len = max(max_len, length(steps))
        
        run_avg_flux = cumsum(flux_samples) ./ steps
        res_flux = abs.(run_avg_flux .- exact_flux) .+ 1e-20 # למניעת log(0)

        # --- עיבוד Correlation ---
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

    # --- הוספת קו תיאורטי 1/sqrt(N) ---
    if max_len > 0
        xs = 100:max_len 
        theory = 0.5 ./ sqrt.(xs) 
        plot!(p2, xs, theory, label="~ 1/√N", color=:black, linestyle=:dot, lw=2)
        plot!(p4, xs, theory, label="~ 1/√N", color=:black, linestyle=:dot, lw=2)
    end

    l = @layout [a b; c d]
    p_final = plot(p1, p2, p3, p4, layout=l, size=(1200, 900), margin=8Plots.mm)
    
    out_name = @sprintf("%s_beta_%.1f_h_%.2f.png", out_prefix, beta, h)
    savefig(p_final, joinpath(OUTPUT_DIR, out_name))
end

# פונקציה להשוואה ישירה: Sequential vs Vectorized (עם קו ED)
function plot_compare_seq_vec(beta, h, eps_inv, path_seq, path_vec, exact_flux, exact_corr)
    println("  Comparing Seq vs Vec vs ED (eps=1/$eps_inv)...")
    
    d_seq = load(path_seq)["measurements"]
    d_vec = load(path_vec)["measurements"]
    
    # חישוב ממוצעים רצים
    f_seq = cumsum(d_seq.avg_spatial_flux) ./ (1:length(d_seq.avg_spatial_flux))
    f_vec = cumsum(d_vec.avg_spatial_flux) ./ (1:length(d_vec.avg_spatial_flux))
    
    c_seq = cumsum(d_seq.flux_corr) ./ (1:length(d_seq.flux_corr))
    c_vec = cumsum(d_vec.flux_corr) ./ (1:length(d_vec.flux_corr))
    
    # גרף Flux
    p1 = plot(title="Flux Comparison (β=$beta, ε=1/$eps_inv)", ylabel="<Flux>", xlabel="Steps")
    plot!(p1, f_seq, label="Sequential", lw=2, alpha=0.8, color=:blue)
    plot!(p1, f_vec, label="Vectorized", lw=2, alpha=0.8, linestyle=:dash, color=:orange)
    hline!(p1, [exact_flux], label="ED Exact", color=:black, linestyle=:dot, lw=2)
    
    # גרף Correlation
    p2 = plot(title="Corr Comparison (β=$beta, ε=1/$eps_inv)", ylabel="<Corr>", xlabel="Steps")
    plot!(p2, c_seq, label="Sequential", lw=2, alpha=0.8, color=:blue)
    plot!(p2, c_vec, label="Vectorized", lw=2, alpha=0.8, linestyle=:dash, color=:orange)
    hline!(p2, [exact_corr], label="ED Exact", color=:black, linestyle=:dot, lw=2)
    
    p_final = plot(p1, p2, layout=(2,1), size=(800, 800), margin=5Plots.mm)
    
    out_name = @sprintf("compare_sv_beta_%.1f_eps_%d_h_%.2f.png", beta, eps_inv, h)
    savefig(p_final, joinpath(OUTPUT_DIR, out_name))
end

# ==========================================
# 4. Main Loop
# ==========================================

function main()
    println("Scanning simulation directory: $SIM_DIR")
    if !isdir(SIM_DIR)
        error("Simulation dir not found!")
    end

    # --- מיפוי הקבצים ---
    registry_seq = []
    registry_vec = []

    for f in readdir(SIM_DIR)
        if endswith(f, ".jld2")
            parsed = parse_sim_filename(f)
            if parsed !== nothing
                path = joinpath(SIM_DIR, f)
                beta, eps_inv, h, is_vec = parsed
                
                if is_vec
                    push!(registry_vec, (path, beta, eps_inv, h))
                else
                    push!(registry_seq, (path, beta, eps_inv, h))
                end
            end
        end
    end
    println("Found $(length(registry_seq)) sequential files.")
    println("Found $(length(registry_vec)) vectorized files.")

    # --- טעינת ED ---
    println("Loading ED data from: $ED_FILE")
    if !isfile(ED_FILE)
        error("ED file not found!")
    end
    ed_data = load(ED_FILE)
    hs_vector = ed_data["hs"]

    # --- לולאת עיבוד ---
    for k in keys(ed_data)
        if k == "hs" continue end
        beta_val = tryparse(Float64, k)
        if beta_val === nothing continue end

        beta_group = ed_data[k]
        if !haskey(beta_group, "F_avg") || !haskey(beta_group, "F_corr") continue end
        
        f_avgs = beta_group["F_avg"]
        f_corrs = beta_group["F_corr"]

        for (i, h_val) in enumerate(hs_vector)
            exact_flux = f_avgs[i]
            exact_corr = f_corrs[i]
            
            # סינון קבצים רלוונטיים
            rel_seq = [(r[1], r[3]) for r in registry_seq if isapprox(r[2], beta_val) && isapprox(r[4], h_val)]
            rel_vec = [(r[1], r[3]) for r in registry_vec if isapprox(r[2], beta_val) && isapprox(r[4], h_val)]

            # 1. גרפים לקבצים רגילים (Sequential) מול ED
            plot_group_vs_ed(beta_val, h_val, exact_flux, exact_corr, rel_seq, "seq")

            # 2. גרפים לקבצים וקטוריים (Vectorized) מול ED
            plot_group_vs_ed(beta_val, h_val, exact_flux, exact_corr, rel_vec, "vec")
            
            # 3. גרף השוואה ישירה (Sequential vs Vectorized)
            for (path_s, eps_s) in rel_seq
                for (path_v, eps_v) in rel_vec
                    if eps_s == eps_v
                        plot_compare_seq_vec(beta_val, h_val, eps_s, path_s, path_v, exact_flux, exact_corr)
                    end
                end
            end
        end
    end
    
    println("\nDone! All plots saved to $OUTPUT_DIR")
end

main()