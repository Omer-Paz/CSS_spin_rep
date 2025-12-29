using GraphRecipes
using Plots
using Graphs
using Random
using SparseArrays
using LinearAlgebra
using Arpack
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
    # We find the sequence [min_val, neighbor1, other, neighbor2]
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

function build_H_in_sector(flux_list,star_list,ne,hx,Kz)
    dim = 0
    states = []
    for i=0:(2^ne - 1)
        state = digits(i,base=2,pad=ne)
        sector = true
        for star in star_list
            star_val = prod([-2*(state[j]-0.5) for j in star])
            if star_val == -1
                sector = false
            end
        end
        if sector == true
            dim += 1
            state_vec = -2 .* (state .- 0.5)
            push!(states,state_vec)
        end
    end
    H = spzeros(dim, dim)
    for (i,state) in enumerate(states)
        diag_energy = -hx * sum(state)
        H[i,i] = diag_energy
        for plaq in flux_list
            new_state = copy(state)
            for p in plaq
                new_state[p] *= -1
            end
            new_state_ind = findfirst(isequal(new_state),states)
            H[i,new_state_ind] += -Kz
        end
    end
    return states,H
end

function flux_op(states,fluxes,N_F)
    flux_op = zeros(length(states),length(states))
    for (i,state) in enumerate(states)
        for plaq in fluxes
            new_state = copy(state)
            for p in plaq
                new_state[p] *= -1
            end
            new_state_ind = findfirst(isequal(new_state),states)
            flux_op[i,new_state_ind] += 1
        end
    end
    return flux_op ./ N_F
end

function flux_corr_op(states,fluxes,l,m)
    flux_op = zeros(length(states),length(states))
    for (i,state) in enumerate(states)
        flux_prod = vcat(fluxes[l], fluxes[m])
        new_state = copy(state)
        for p in flux_prod
            new_state[p] *= -1
        end
        new_state_ind = findfirst(isequal(new_state),states)
        flux_op[i,new_state_ind] += 1
    end
    return flux_op
end

function calc_thermal_avg(evals, evecs, Op, beta)
    if beta == Inf
        psi_gs = evecs[:, 1]
        psi_gs = psi_gs / norm(psi_gs)
        return dot(psi_gs, Op, psi_gs)
    else
        E_min = minimum(evals)
        weights = exp.(-beta .* (evals .- E_min))
        Z = sum(weights)
        expectations = [evecs[:, i]' * Op * evecs[:, i] for i in 1:length(evals)]
        return sum(expectations .* weights) / Z
    end
end

function build_graph_data(g)
    four_cycles = find_all_4_cycles(g)
    flux_list = [[(min(src(e),dst(e)),max(src(e),dst(e))) for e in f] for f in four_cycles]
    flux_list=[sort(f, by=first) for f in flux_list]
    flux_list=[sort(f, by=x->x[2]) for f in flux_list]
    flux_list=unique(flux_list)
    edge_list = [(min(src(e),dst(e)),max(src(e),dst(e))) for e in edges(g)]
    fluxes = []
    for flux in flux_list
        f = []
        for e in flux
            ind = findfirst(isequal(e),edge_list)
            push!(f,ind)
        end
        push!(fluxes,f)
    end
    star_list = [Int[] for _ in 1:nv(g)]
    for (i, edge_tuple) in enumerate(edge_list)
        u, v = edge_tuple
        push!(star_list[u], i)
        push!(star_list[v], i)
    end
    return flux_list,fluxes,star_list
end

function run_ED(g,K_z,hs,betas,to_save,folder_name,l=1,m=6)
    flux_list,fluxes,star_list = build_graph_data(g)
    N_F = length(fluxes)
    ED_data = Dict()
    ED_data["hs"] = hs
    for beta in betas 
        ED_data["$beta"] = Dict("F_avg"=>[], "F_corr"=>[],"E"=>[])
    end
    for h in hs
        states,H = build_H_in_sector(fluxes,star_list,ne(g),h,K_z)
        f_op = flux_op(states,fluxes,N_F)
        f_corr_op = flux_corr_op(states,fluxes,l,m)    
        vals, vecs = eigs(H; nev=length(states), which=:SR)
        for beta in betas
            avg_flux = calc_thermal_avg(vals,vecs,f_op,beta)
            flux_corr = calc_thermal_avg(vals,vecs,f_corr_op,beta)
            H_per_flux = H ./ N_F
            energy = calc_thermal_avg(vals,vecs,H_per_flux,beta)
            push!(ED_data["$beta"]["F_avg"],avg_flux)
            push!(ED_data["$beta"]["F_corr"],flux_corr)
            push!(ED_data["$beta"]["E"],energy)
        end
    end
    if to_save
        save("/Users/omerp/Cayley_CSS-1/ED/$folder_name/$folder_name"*"_ED.jld2",ED_data)
    end
end
# ----------------------------------------------------------------------------------------------------

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
#@show ne(g)
#Random.seed!(1234)
#graphplot(g, method=:spring, curves=false,markersize=0.1)

# ----------------------------------------------------------------------------------------------------


K_z = 1.0
hs = range(0.0,1.0,11)
Fs = []
corrs = []
betas = [0.5,1.0,2.0,4.0,8.0,16.0]
to_save = true
folder_name = basename(@__DIR__)
run_ED(g,K_z,hs,betas,to_save,folder_name)

#betas = [0.5,1.0,2.0,4.0,8.0,16.0]
#fig = plot()
#for beta in betas
#    hs = d["hs"]
#    F_avg = d["$beta"]["F_avg"]
#    plot!(hs,F_avg,label="beta=$beta")
#end
#display(fig)
