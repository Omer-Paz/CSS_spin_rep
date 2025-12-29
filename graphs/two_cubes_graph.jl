using GraphRecipes
using Plots
using Graphs
using Random
using JLD2

function canonical_cycle_rep(cycle::Vector{Int})
    if length(cycle) != 4
        error("Cycle must have 4 vertices")
    end

    # Find the index of the minimum vertex
    min_val = cycle[1]
    min_idx = 1
    for i in 2:4
        if cycle[i] < min_val
            min_val = cycle[i]
            min_idx = i
        end
    end

    # Shift the cycle so it starts with the minimum vertex.
    local shifted_cycle::Vector{Int}
    if min_idx == 1
        shifted_cycle = cycle
    elseif min_idx == 2
        shifted_cycle = [cycle[2], cycle[3], cycle[4], cycle[1]]
    elseif min_idx == 3
        shifted_cycle = [cycle[3], cycle[4], cycle[1], cycle[2]]
    else # min_idx == 4
        shifted_cycle = [cycle[4], cycle[1], cycle[2], cycle[3]]
    end

    if shifted_cycle[2] < shifted_cycle[4]
        return shifted_cycle
    else
        return [shifted_cycle[1], shifted_cycle[4], shifted_cycle[3], shifted_cycle[2]]
    end
end

function find_all_4_cycles(g::SimpleGraph)
    cycles_set = Set{Set{Edge}}()

    for u in vertices(g)
        for v in neighbors(g, u)
            for w in neighbors(g, v)
                if w == u
                    continue
                end

                for x in neighbors(g, w)
                    if x == v || x == u
                        continue
                    end

                    if has_edge(g, x, u)
                        edge_set = Set([
                            Edge(u, v),
                            Edge(v, w),
                            Edge(w, x),
                            Edge(x, u)
                        ])
                        push!(cycles_set, edge_set)
                    end
                end
            end
        end
    end

    return cycles_set
end

# --- פונקציות עזר לייצור מיפויים ---

function generate_E_τ(edge_list, flux_list)
    E_τ = Vector{Vector{Int}}(undef, length(edge_list))
    for (n, edge) in enumerate(edge_list)
        E_τ[n] = [i for (i,f) in enumerate(flux_list) if edge in f]
    end
    return E_τ
end

function generate_flux_to_edge_indices(flux_list_nodes, edge_list)
    edge_lookup = Dict{Tuple{Int,Int}, Int}()
    for (i, e) in enumerate(edge_list)
        edge_lookup[e] = i
    end

    flux_to_edge = Vector{Vector{Int}}(undef, length(flux_list_nodes))
    for (i, flux) in enumerate(flux_list_nodes)
        flux_to_edge[i] = [edge_lookup[e] for e in flux]
    end
    return flux_to_edge
end

# הפונקציה החדשה שביקשת
function generate_vert_to_edge(g::SimpleGraph, edge_list)
    # אתחול רשימה של רשימות, אחת לכל צומת
    v_to_e = [Int[] for _ in 1:nv(g)]
    
    # מעבר על רשימת הקשתות ושיוך האינדקס לצמתים המתאימים
    for (i, (u, v)) in enumerate(edge_list)
        push!(v_to_e[u], i)
        push!(v_to_e[v], i)
    end
    return v_to_e
end


# --- בניית הגרף ---

g = SimpleGraph(12)
add_edge!(g,1,2)
add_edge!(g,2,3)
add_edge!(g,3,4)
add_edge!(g,4,1)

add_edge!(g,1+4,2+4)
add_edge!(g,2+4,3+4)
add_edge!(g,3+4,4+4)
add_edge!(g,4+4,1+4)

add_edge!(g,1+8,2+8)
add_edge!(g,2+8,3+8)
add_edge!(g,3+8,4+8)
add_edge!(g,4+8,1+8)

add_edge!(g,1,5)
add_edge!(g,5,9)
add_edge!(g,2,6)
add_edge!(g,6,10)
add_edge!(g,3,7)
add_edge!(g,7,11)
add_edge!(g,4,8)
add_edge!(g,8,12)

# --- עיבוד הנתונים ---

# יצירת רשימת פלאקטים
four_cycles = find_all_4_cycles(g)
flux_list = [[(min(src(e),dst(e)),max(src(e),dst(e))) for e in f] for f in four_cycles]
flux_list=[sort(f, by=first) for f in flux_list]
flux_list=[sort(f, by=x->x[2]) for f in flux_list]
flux_list=unique(flux_list)

# יצירת רשימת קשתות
edge_list = [(min(src(e),dst(e)),max(src(e),dst(e))) for e in edges(g)]

# יצירת המיפויים
E_τ = generate_E_τ(edge_list, flux_list)
flux_to_edge_indices = generate_flux_to_edge_indices(flux_list, edge_list)
vert_to_edge = generate_vert_to_edge(g, edge_list) # <--- קריאה לפונקציה החדשה

# --- שמירה ---

d = Dict()
d["N_vertices"] = nv(g)          # הוספתי גם את זה, שימושי
d["N_edges"]    = length(edge_list)
d["N_fluxes"]   = length(flux_list)

d["edge_to_vert"] = edge_list
d["edge_to_flux"] = E_τ
d["flux_to_edge"] = flux_to_edge_indices
d["vert_to_edge"] = vert_to_edge # <--- הוספה למילון

save_path = @__DIR__
save(joinpath(save_path, "two_cubes.jld2"), d)

println("Saved geometry including vert_to_edge to two_cubes.jld2")