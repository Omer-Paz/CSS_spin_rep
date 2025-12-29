using Nemo
using Graphs, DataStructures, Plots, GraphRecipes, Random
using NetworkLayout
using JLD2

function canonical_rep(M)
    F = parent(M[1,1])
    for val in M
        if val != 0
            return 1/val .* M
        end
    end
end

function general_linear(n)
    F,a = finite_field(2,2*n,"a")
    GL = []
    for x in F, y in F, z in F, t in F
        det = x*t-y*z
        if det==0
            continue
        end
        push!(GL,[x y; z t])
    end
    return GL
end

function special_linear(n)
    GL = general_linear(n)
    for i=1:length(GL)
        M = canonical_rep(GL[i])
        GL[i] = M
    end
    return unique(GL)
end

function morgerstern_gen(n)
    F,a = finite_field(2,2*n,"a")
    i = a^2
    for x in F
        if x^2 + x == 1
            i = x
        end
    end
    s_0 = canonical_rep([a^0 i; a*(1+i) a^0])
    s_1 = canonical_rep([a^0 a^0+i; a*i a^0])
    s_2 = canonical_rep([a^0 a^0; a a^0])
    return [s_0,s_1,s_2], [canonical_rep(s_0*s_1),canonical_rep(s_0*s_2),canonical_rep(s_1*s_0),canonical_rep(s_2*s_0)]
end


function cayley_graph(G,S,lr)
    index_map = Dict(M => i for (i, M) in enumerate(G))  # Map elements to indices
    g = SimpleGraph(length(G))  # Create empty graph
    for M in G
        u = index_map[M]
        for gen in S
            if lr == 1
                N = canonical_rep(M * gen)  # Matrix multiplication in PSL(2, q)
            elseif lr == -1
                N = canonical_rep(gen * M)
            end
            v = index_map[N]  # Get vertex index
            add_edge!(g, u, v)
        end
    end

    return g
end

function lr_cayley_complex(G,A,B)
    flux_list = []
    index_map = Dict(M => i for (i, M) in enumerate(G))  # Map elements to indices
    g = SimpleGraph(length(G))  # Create empty graph
    for M in G
        u = index_map[M]
        for a∈A
            N = canonical_rep(a * M)  # Matrix multiplication in PSL(2, q)
            v = index_map[N]  # Get vertex index
            add_edge!(g, u, v)
        end
        for b∈B
            N = canonical_rep(M * b)
            v=index_map[N]
            add_edge!(g, u, v)
        end
        for a∈A, b∈B 
            f = [] 
            M_1 = canonical_rep(a*M) 
            M_2 = canonical_rep(M*b) 
            M_3 = canonical_rep(a*M*b) 
            u_1 = index_map[M_1] 
            u_2 = index_map[M_2] 
            u_3 = index_map[M_3] 
            f = [sort([u,u_1]),sort([u_1,u_3]),sort([u_3,u_2]),sort([u_2,u])] 
            f = sort(f, by=x -> x[2])
            f = sort(f, by=x -> x[1])
            f = [(f[1][1],f[1][2]),(f[2][1],f[2][2]),(f[3][1],f[3][2]),(f[4][1],f[4][2])]
            push!(flux_list,f)
        end
    end
    flux_list = unique(flux_list)
    flux_list = sort(flux_list, by=x -> x[4])
    flux_list = sort(flux_list, by=x -> x[3])
    flux_list = sort(flux_list, by=x -> x[2])
    flux_list = sort(flux_list, by=x -> x[1])
    edge_list = []
    for e in collect(edges(g))
        v_1 = src(e)
        v_2 = dst(e)
        if v_1 > v_2
            push!(edge_list,(v_2,v_1))
        else
            push!(edge_list,(v_1,v_2))
        end
    end
    return g,flux_list,edge_list
end

#=function generate_E_τ(edge_list,flux_list::Vector{Tuple{Int64,Int64}}) # for every edge e, find the fluxes (cubes centered above fluxes) containing e
    N_F_τ = length(edge_list)
    E_τ = []
    for n=1:N_F_τ # go over edges
        l = []
        edge = edge_list[n]
        # if edge in flux_list[k], push k to l
        for f in flux_list
            for f_edge in f
                if f_edge==edge
                    println("typeof F=",typeof(f))
                    println("typeof F list = ",typeof(flux_list))
                    f_index = findfirst(f,flux_list)
                    push!(l,f_index)
                end
            end
        end
        push!(E_τ,l)
    end
    return E_τ
end=#

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

#=function visualize_graph(g)
    Random.seed!(1234)
    graphplot(g, method=:spring, curves=false,markersize=0.1)
end

function visualize_two_graphs(g1, g2; line_width=1,a=1, color1=:red, color2=:black , stroke_width=2,mk_size=5)
    Random.seed!(1234)
    layout = NetworkLayout.spring(g1)  # Vertices location
    xs = [p[1] for p in layout]
    ys = [p[2] for p in layout]

    # Plot first graph’s edges
    plt = plot(legend=false, axis=false, size=(700,700))
    for e in edges(g1)
        u, v = src(e), dst(e)
        plot!([xs[u], xs[v]], [ys[u], ys[v]], color=color1, alpha=a, lw=line_width)
    end

    # Overlay second graph’s edges
    for e in edges(g2)
        u, v = src(e), dst(e)
        plot!([xs[u], xs[v]], [ys[u], ys[v]], color=color2, alpha=a, lw=line_width)
    end

    # Add vertices
    scatter!(xs, ys, color=:black, markersize=mk_size, markercolor=:blue,markerstrokewidth=stroke_width,grid=false)

    return plt
end

function get_check_matrices(g, z_checks)
    num_qubits = ne(g) 
    num_z_checks = length(z_checks)
    H_z = zeros(Int, num_z_checks, num_qubits)
    edge_map = Dict{Tuple{Int, Int}, Int}()
    for (idx, edge) in enumerate(edges(g))
        u, v = edge.src, edge.dst
        edge_map[(u, v)] = idx
        edge_map[(v, u)] = idx
    end
    for (i, check) in enumerate(z_checks) # 'i' is the row index (check index)
        for edge_tuple in check
            if haskey(edge_map, edge_tuple)
                qubit_idx = edge_map[edge_tuple] # 'qubit_idx' is the col index
                H_z[i, qubit_idx] = 1
            else
                # This warning is helpful for debugging
                println("Warning: Edge $edge_tuple in Z check $i not found in graph.")
            end
        end
    end
    H_x = Matrix{Int}(incidence_matrix(g))
    return H_z,H_x
end



function kernel_dim(A)
    A = copy(A) .% 2  # work mod 2
    m, n = size(A)
    rank = 0
    row = 1
    for col in 1:n
        # find pivot
        pivot = findfirst(r -> A[r, col] == 1, row:m)
        if pivot !== nothing
            pivot += row - 1
            # swap pivot row up
            if pivot != row
                A[row, :], A[pivot, :] = A[pivot, :], A[row, :]
            end
            # eliminate below and above
            for r in 1:m
                if r != row && A[r, col] == 1
                    A[r, :] = (A[r, :] .⊻ A[row, :])  # XOR = addition mod 2
                end
            end
            rank += 1
            row += 1
            if row > m
                break
            end
        end
    end

    # kernel dimension = number of columns - rank
    return n - rank
end
=#
n=1
SL = special_linear(n);
println("Size SL = ",size(SL))
A,B = morgerstern_gen(n);
r_cayley = cayley_graph(SL,A,1);
l_cayley = cayley_graph(SL,B,-1);
lr_cayley,flux_list,edge_list = lr_cayley_complex(SL,A,B);
E_τ = generate_E_τ(edge_list, flux_list)
flux_to_edge_indices = generate_flux_to_edge_indices(flux_list, edge_list)
vert_to_edge = generate_vert_to_edge(lr_cayley, edge_list) 

d = Dict()
d["N_vertices"] = nv(lr_cayley)         
d["N_edges"]    = length(edge_list)
d["N_fluxes"]   = length(flux_list)

d["edge_to_vert"] = edge_list
d["edge_to_flux"] = E_τ
d["flux_to_edge"] = flux_to_edge_indices
d["vert_to_edge"] = vert_to_edge 
d["l_cayley"] = l_cayley
d["r_cayley"] = r_cayley
d["lr_cayley"] = lr_cayley
d["A"] = A
d["B"] = B
field_size = 2^(2*n)
file_path = joinpath(@__DIR__,"psl_2_$(field_size)_graph.jld2")
save(file_path,d)
d

#H_z, H_x = get_check_matrices(lr_cayley,flux_list)
#k_x = kernel_dim(H_x)
#k_z = kernel_dim(H_z)
#println(k_x)
#println(k_z)
#println("nm edges = ",length(edge_list))
#E_τ = generate_E_τ(edge_list,flux_list);
#C_to_E = generate_C_to_E(E_τ,flux_list);
#d = Dict()
#field_size = 2^(2*n)
#save("psl_2_$field_size"*"_hg.jld2",d)
#d
#println(flux_list)
#println(dst(collect(edges(lr_cayley))[1]))
#visualize_graph(lr_cayley)
#visualize_graph(l_cayley)
#visualize_two_graphs(l_cayley,r_cayley)





























#println(A[1] in SL)
#println(A[2] in SL)
#println(A[3] in SL)
#println(B[1] in SL)
#println(B[2] in SL)
#println(B[3] in SL)
#println(B[4] in SL)
#println(canonical_rep(A[1]^2))
#println(canonical_rep(A[2]^2))
#println(canonical_rep(A[3]^2))
#println(canonical_rep(B[1]^2))
#println(canonical_rep(B[2]^2))
#println(canonical_rep(B[3]^2))
#println(canonical_rep(B[4]^2))