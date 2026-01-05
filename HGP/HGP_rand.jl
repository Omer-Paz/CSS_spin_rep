using Random
using SparseArrays
using LinearAlgebra
using JLD2

function rank_gf2(H::SparseMatrixCSC)
    m, n = size(H)
    A = Matrix{Int}(H) .% 2 
    r = 0
    pivot_row = 1
    for j in 1:n
        pivot_row > m && break
        idx = findfirst(x -> x == 1, A[pivot_row:m, j])
        if idx !== nothing
            idx += pivot_row - 1
            A[pivot_row, :], A[idx, :] = A[idx, :], A[pivot_row, :]
            for i in 1:m
                if i != pivot_row && A[i, j] == 1
                    A[i, :] .= (A[i, :] .+ A[pivot_row, :]) .% 2
                end
            end
            r += 1
            pivot_row += 1
        end
    end
    return r
end

function classical_distance(H::SparseMatrixCSC)
    m, n = size(H)
    min_dist = n
    for i in 1:(2^n - 1)
        x = digits(i, base=2, pad=n)
        if iszero((H * x) .% 2)
            w = sum(x)
            if 0 < w < min_dist
                min_dist = w
            end
        end
    end
    return min_dist
end

function HPG_analysis(H1, H2, HX, HZ)
    n1, n2 = size(H1, 2), size(H2, 2)
    m1, m2 = size(H1, 1), size(H2, 1)
    N = n1*n2 + m1*m2
    
    rkX = rank_gf2(HX)
    rkZ = rank_gf2(HZ)
    k = N - rkX - rkZ
    
    d1 = classical_distance(H1)
    d2 = classical_distance(H2)
    d1T = classical_distance(sparse(H1'))
    d2T = classical_distance(sparse(H2'))
    d_lower_bound = min(d1, d2, d1T, d2T)
    
    return k, d_lower_bound
end

function gallager_ldpc_construction(n::Int, dv::Int, dc::Int)
    @assert isinteger(n * dv ÷ dc)
    m = n * dv ÷ dc  # Number of check nodes
    bits_inds = repeat(1:n, inner=dv)
    check_inds = repeat(1:m, outer=dc)
    @assert length(bits_inds) == length(check_inds)
    shuffle!(check_inds)
    edges = collect(zip(bits_inds, check_inds))
    counts = Dict{Tuple{Int, Int},Int}()
    for e in edges
        counts[e] = get(counts,e,0) + 1
    end
    doubled_edges = []
    for e in edges
        if counts[e] > 1
            push!(doubled_edges, e)
        end
    end
    if !isempty(doubled_edges)
        return gallager_ldpc_construction(n, dv, dc)
    end
    adjecancy = spzeros(Bool, m, n)
    for (b, c) in edges
        adjecancy[c, b] = true
    end
    return adjecancy
end

function HPG_product(H_1::SparseMatrixCSC,H_2::SparseMatrixCSC)
    m1, n1 = size(H_1)
    m2, n2 = size(H_2)
    In1 = sparse(I, n1, n1)
    Im1 = sparse(I, m1, m1)
    Im2 = sparse(I, m2, m2)
    In2 = sparse(I, n2, n2)
    H_X = hcat(kron(H_1, In2), kron(Im1, H_2'))
    H_Z = hcat(kron(In1, H_2), kron(H_1', Im2))
    CSS_condition = (H_Z * H_X') .% 2
    @assert iszero(CSS_condition)
    return H_X, H_Z 
end

function create_flux_to_edge(H_Z)
    num_z_checks = size(H_Z, 1)
    return [findall(!iszero, H_Z[i, :]) for i in 1:num_z_checks]
end

function create_edge_to_flux(H_Z)
    num_qubits = size(H_Z,2)
    return [findall(!iszero, H_Z[:, j]) for j in 1:num_qubits]
end

function create_vert_to_edge(H_X)
    num_x_checks = size(H_X, 1)
    return [findall(!iszero, H_X[i, :]) for i in 1:num_x_checks]
end

function create_edge_to_vert(H_X)
    num_qubits = size(H_X,2)
    return [findall(!iszero, H_X[:, j]) for j in 1:num_qubits]
end


n_bits = 10
dv = 3
dc = 6
H_1 = gallager_ldpc_construction(n_bits, dv, dc)
H_2 = gallager_ldpc_construction(n_bits, dv, dc)
H_X, H_Z = HPG_product(H_1, H_2)
k,d_lower_bound =  HPG_analysis(H_1, H_2, H_X, H_Z)


vertices = [findall(x->x!=0, H_X[:,j]) for j in 1:size(H_X, 2)] # X-checks
edges = collect(1:size(H_X, 2)) # spatial edges

flux_to_edge =  create_flux_to_edge(H_Z) # fluxes over spatial edges
edge_to_flux = create_edge_to_flux(H_Z)
edge_to_vert = create_edge_to_vert(H_X)
vert_to_edge = create_vert_to_edge(H_X)

d = Dict()
d["N_edges"] = size(H_X, 2)
d["N_vertices"] = size(H_X, 1) # vertex == check
d["N_fluxes"] = size(H_Z, 1)
d["flux_to_edge"] = flux_to_edge
d["edge_to_flux"] = edge_to_flux
d["vert_to_edge"] = vert_to_edge
d["edge_to_vert"] = edge_to_vert
d["k"] = k
d["d_lower_bound"] = d_lower_bound
save_path = joinpath(@__DIR__,"HGP_test.jld2")
save(save_path,d)