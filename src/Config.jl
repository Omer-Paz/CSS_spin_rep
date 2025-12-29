mutable struct SimConfig
    # == State (Dynamic) ==
    ising_config_G::Matrix{Int8}    # Size: (N_edges, M)
    ising_config_tau::Matrix{Int8}  # Size: (N_vertices, M)
    F_G::Matrix{Int8}     # Spatial Fluxes (Plaquettes). Size: (N_fluxes, M)
    F_tau::Matrix{Int8}   # Temporal Fluxes. Size: (N_edges, M) 
    # == Geometry (Static - Loaded from JLD2) ==
    edge_to_vert::Vector{Tuple{Int, Int}}  
    edge_to_flux::Vector{Vector{Int}}      
    flux_to_edge::Vector{Vector{Int}}     
    vert_to_edge::Vector{Vector{Int}} 
    # == Buffers ==
    temp_dEs::Vector{Float64}      
    temp_randoms::Vector{Float64}  
    G_dEs::Vector{Float64}
    G_randoms::Vector{Float64}
    # == Constructor ==
    function SimConfig(geometry_dict::Dict, M::Int)
        N_edges = geometry_dict["N_edges"]
        N_vertices = geometry_dict["N_vertices"] 
        N_fluxes = geometry_dict["N_fluxes"]
        config_G = ones(Int8, N_edges, M)
        config_tau = ones(Int8, N_vertices, M)
        flux_G = ones(Int8, N_fluxes, M)
        flux_tau = ones(Int8, N_edges, M)
        temp_dEs = zeros(Float64, M)
        temp_randoms = zeros(Float64, M)
        G_dEs = zeros(Float64, M)
        G_randoms = zeros(Float64, M)
        new(config_G, config_tau, flux_G, flux_tau,
            geometry_dict["edge_to_vert"],
            geometry_dict["edge_to_flux"],
            geometry_dict["flux_to_edge"],
            geometry_dict["vert_to_edge"],
            temp_dEs,
            temp_randoms,
            G_dEs,
            G_randoms)
    end
end