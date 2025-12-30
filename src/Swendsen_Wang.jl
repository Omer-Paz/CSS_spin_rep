using Random
using Graphs
struct TannerGraph
    spatial_spin_to_checks::Dict{Tuple{Int, Int}, Vector{Tuple{Symbol, Int, Int}}}
    temporal_spin_to_checks::Dict{Tuple{Int, Int}, Vector{Tuple{Symbol, Int, Int}}}
    check_to_spins::Dict{Tuple{Symbol, Int, Int}, Vector{Tuple{Symbol, Int, Int}}}
    spin_degrees::Dict{Tuple{Symbol, Int, Int}, Int} # For peeling algorithm
end

function TannerGraph()
    return TannerGraph(
        Dict{Tuple{Int, Int}, Vector{Tuple{Symbol, Int, Int}}}(),
        Dict{Tuple{Int, Int}, Vector{Tuple{Symbol, Int, Int}}}(),
        Dict{Tuple{Symbol, Int, Int}, Vector{Tuple{Symbol, Int, Int}}}(),
        Dict{Tuple{Symbol, Int, Int}, Int}()
    )
end

function form_clusters(conf::SimConfig, params::SimParams)
    p_G = params.p_G
    p_tau = params.p_tau
    active_G_checks = Tuple{Int, Int}[]
    active_tau_checks = Tuple{Int, Int}[]
    
    M = params.M
    N_fluxes = size(conf.F_G, 1)
    N_edges = size(conf.F_tau, 1)
    for m in 1:M
        for p in 1:N_fluxes
            if conf.F_G[p, m] == 1 && rand(params.rng) < p_G
                push!(active_G_checks, (p, m))
            end
        end
    end
    for m in 1:M
        for s in 1:N_edges
            if conf.F_tau[s, m] == 1 && rand(params.rng) < p_tau
                push!(active_tau_checks, (s, m))
            end
        end
    end

    return active_G_checks, active_tau_checks
end


function add_connection!(graph::TannerGraph,check_key, spin_type, spin_idx, spin_time)
    spin_key = (spin_type, spin_idx, spin_time)
    
    if !haskey(graph.check_to_spins, check_key)
        graph.check_to_spins[check_key] = []
    end
    push!(graph.check_to_spins[check_key], spin_key)
    
    # 2. הוספת הבדיקה לרשימת הבדיקות של הספין + עדכון דרגה
    if spin_type == :spatial
        dict = graph.spatial_spin_to_checks
        raw_key = (spin_idx, spin_time)
    else
        dict = graph.temporal_spin_to_checks
        raw_key = (spin_idx, spin_time)
    end
    if !haskey(dict, raw_key)
        dict[raw_key] = []
    end
    push!(dict[raw_key], check_key)
    graph.spin_degrees[spin_key] = get(graph.spin_degrees, spin_key, 0) + 1
end

    
function build_tanner_graph(conf::SimConfig, M::Int, active_G_checks::Vector{Tuple{Int, Int}}, active_tau_checks::Vector{Tuple{Int, Int}})
    graph = TannerGraph()

    # --- הוספת פלאקטים מרחביים ---
    for (p_idx, m) in active_G_checks
        check_key = (:G, p_idx, m)
        # שימוש בגיאומטריה הקיימת: flux_to_edge
        spatial_edges = conf.flux_to_edge[p_idx]

        for edge_idx in spatial_edges
            # קשתות מרחביות נמצאות באותו זמן m כמו הפלאקט
            add_connection!(graph,check_key, :spatial, edge_idx, m)
        end
    end
    for (s_idx, m) in active_tau_checks
        check_key = (:tau, s_idx, m)
        m_next = (m == M) ? 1 : m + 1
        add_connection!(graph,check_key, :spatial, s_idx, m)
        add_connection!(graph,check_key, :spatial, s_idx, m_next)

        # 2. קשתות זמניות (temporal edges / vertices)
        # אלו הוורטיקלים שמחוברים לקצוות הקשת s
        # משתמשים ב-edge_to_vert מהקונפיג
        u, v = conf.edge_to_vert[s_idx]

        # קשתות זמניות נמצאות "בין" השכבות, נהוג לשייך אותן לאינדקס m
        add_connection!(graph,check_key, :temporal, u, m)
        add_connection!(graph,check_key, :temporal, v, m)
    end
    return graph
end


# ==========================================
# חלק 3: Peeling Decoder
# ==========================================

"""
מבצע את אלגוריתם ה-Peeling בשני מעברים (Forward & Backward)
כדי לדגום קונפיגורציה חוקית מגרעין האילוצים.
"""
function run_peeling_decoder!(conf::SimConfig, params::SimParams, graph::TannerGraph)
    # 1. אתחול אקראי של כל המערכת
    # זה קובע ערכים למשתנים החופשיים ומנחש ערכים ראשוניים לתלויים
    rand!(params.rng, conf.ising_config_G, [-1, 1])
    rand!(params.rng, conf.ising_config_tau, [-1, 1])
    
    # 2. Forward Pass - בניית מחסנית התיקונים
    # מחזיר רשימה של (spin_type, index, time, check_to_satisfy)
    correction_stack = forward_peeling_pass(graph)
    
    # 3. Backward Pass - תיקון האילוצים
    apply_corrections!(conf, graph, correction_stack)
    
    # עדכון שטפים (Fluxes) שיהיו תואמים לספינים החדשים
    update_all_fluxes!(conf, params)
end

function forward_peeling_pass(graph::TannerGraph)
    # יצירת עותקים ברי-שינוי למעקב אחרי מצב הגרף
    current_degrees = copy(graph.spin_degrees)
    
    # מעקב אחרי בדיקות פעילות (כדי לא למחוק פעמיים)
    active_checks = Set{Tuple{Symbol, Int, Int}}(keys(graph.check_to_spins))
    
    # תור עלים (ספינים עם דרגה 1)
    leaf_queue = Tuple{Symbol, Int, Int}[]
    
    # רשימת כל הספינים הפעילים (כדי לבחור שרירותית במקרה של לולאה)
    active_spins = Set{Tuple{Symbol, Int, Int}}(keys(current_degrees))

    # אתחול התור
    for (spin, deg) in current_degrees
        if deg == 1
            push!(leaf_queue, spin)
        end
    end

    # המחסנית שתחזיק את סדר הפעולות: (spin, check_that_solves_it)
    stack = Vector{Tuple{Tuple{Symbol, Int, Int}, Tuple{Symbol, Int, Int}}}()

    while !isempty(active_checks)
        # מקרה א': יש עלים
        if !isempty(leaf_queue)
            spin = popfirst!(leaf_queue)
            
            # אם הספין כבר טופל (הוסר מהגרף), דלג
            if !haskey(current_degrees, spin)
                continue
            end
            
            # מצא את הבדיקה היחידה שמחוברת לספין הזה
            # (בגלל שהדרגה היא 1, חייבת להיות בדיוק אחת פעילה)
            connected_checks = (spin[1] == :spatial) ? 
                               graph.spatial_spin_to_checks[spin[2:3]] : 
                               graph.temporal_spin_to_checks[spin[2:3]]
            
            target_check = nothing
            for check in connected_checks
                if check in active_checks
                    target_check = check
                    break
                end
            end
            
            # אם לא מצאנו (באג לוגי או שהספין כבר לא רלוונטי), דלג
            if isnothing(target_check)
                delete!(current_degrees, spin)
                delete!(active_spins, spin)
                continue
            end

            # דחוף למחסנית: הספין הזה יפתור את הבדיקה הזו
            push!(stack, (spin, target_check))
            
            # "מחק" את הבדיקה והספין מהגרף
            delete!(active_checks, target_check)
            delete!(current_degrees, spin)
            delete!(active_spins, spin)
            
            # עדכן את השכנים של הבדיקה שנמחקה
            neighbors = graph.check_to_spins[target_check]
            for neighbor_spin in neighbors
                if haskey(current_degrees, neighbor_spin)
                    current_degrees[neighbor_spin] -= 1
                    if current_degrees[neighbor_spin] == 1
                        push!(leaf_queue, neighbor_spin)
                    end
                end
            end

        # מקרה ב': אין עלים (לולאה סגורה) - שבירת סימטריה
        else
            # בחר ספין אקראי שעדיין פעיל
            if isempty(active_spins)
                break # סיימנו
            end
            rand_spin = first(active_spins) # בחירה שרירותית
            
            # מצא בדיקה כלשהי שמחוברת אליו כדי להפוך אותה ל"פתורה" ע"י הספין הזה
            # (אנחנו מכריחים אותו להיות עלה באופן מלאכותי)
            connected_checks = (rand_spin[1] == :spatial) ? 
                               graph.spatial_spin_to_checks[rand_spin[2:3]] : 
                               graph.temporal_spin_to_checks[rand_spin[2:3]]
            
            forced_check = nothing
            for check in connected_checks
                if check in active_checks
                    forced_check = check
                    break
                end
            end
            
            if !isnothing(forced_check)
                # דחוף למחסנית כאילו הוא היה עלה
                push!(stack, (rand_spin, forced_check))
                
                # הסר את הבדיקה (וכך נפתחת הלולאה לאחרים)
                delete!(active_checks, forced_check)
                neighbors = graph.check_to_spins[forced_check]
                for neighbor_spin in neighbors
                    if haskey(current_degrees, neighbor_spin)
                        current_degrees[neighbor_spin] -= 1
                        if current_degrees[neighbor_spin] == 1
                            push!(leaf_queue, neighbor_spin)
                        end
                    end
                end
            end
            
            # הסר את הספין הזה מהמשחק
            delete!(current_degrees, rand_spin)
            delete!(active_spins, rand_spin)
        end
    end
    
    return stack
end

function apply_corrections!(conf::SimConfig, graph::TannerGraph, stack)
    # עוברים על המחסנית מהסוף להתחלה (LIFO)
    while !isempty(stack)
        (spin_key, check_key) = pop!(stack)
        
        # 1. חשב את ה-Parity הנוכחי של הבדיקה
        # (מכפלת הספינים המחוברים אליה)
        neighbors = graph.check_to_spins[check_key]
        current_parity = 1
        
        for (sType, sIdx, sTime) in neighbors
            val = (sType == :spatial) ? 
                  conf.ising_config_G[sIdx, sTime] : 
                  conf.ising_config_tau[sIdx, sTime]
            current_parity *= val
        end
        
        # 2. אם הבדיקה לא מסופקת (מכפלה -1), הפוך את הספין ה"תלוי"
        # הספין הזה נבחר ב-Forward pass להיות האחראי על הבדיקה הזו
        if current_parity == -1
            (sType, sIdx, sTime) = spin_key
            if sType == :spatial
                conf.ising_config_G[sIdx, sTime] *= -1
            else
                conf.ising_config_tau[sIdx, sTime] *= -1
            end
        end
    end
end

# פונקציית עזר לעדכון כל המטריצות F_G ו-F_tau אחרי שינוי הספינים
function update_all_fluxes!(conf::SimConfig, params::SimParams)
    # עדכון שטפים מרחביים
    N_fluxes, M = size(conf.F_G)
    for m in 1:M
        for p in 1:N_fluxes
            edges = conf.flux_to_edge[p]
            prod_val = 1
            for e in edges
                prod_val *= conf.ising_config_G[e, m]
            end
            conf.F_G[p, m] = prod_val
        end
    end
    
    # עדכון שטפים זמניים
    N_edges, M = size(conf.F_tau)
    for m in 1:M
        m_next = (m == M) ? 1 : m + 1
        for s in 1:N_edges
            # פלאקט זמני: s(m) * s(m+1) * tau_u(m) * tau_v(m)
            u, v = conf.edge_to_vert[s]
            
            val = conf.ising_config_G[s, m] * conf.ising_config_G[s, m_next] * conf.ising_config_tau[u, m] * conf.ising_config_tau[v, m]
                  
            conf.F_tau[s, m] = val
        end
    end
end

function run_swendsen_wang_step!(conf::SimConfig, params::SimParams)
    active_G, active_tau = form_clusters(conf, params)
    graph = build_tanner_graph(conf, params.M, active_G, active_tau)
    run_peeling_decoder!(conf, params, graph)
end