using Random
using LoopVectorization

@inline function next_time(m::Int, M::Int)
    return m == M ? 1 : m + 1
end

@inline function prev_time(m::Int, M::Int)
    return m == 1 ? M : m - 1
end

function calc_dE_spatial(conf::SimConfig, params::SimParams, s::Int, m::Int)
    sum_flux_spatial = 0.0
    @inbounds for p_idx in conf.edge_to_flux[s]
        sum_flux_spatial += conf.F_G[p_idx, m]
    end
    
    dE_spatial = 2.0 * params.G_coupling * sum_flux_spatial
    m_prev = prev_time(m, params.M)
    flux_tau_up = conf.F_tau[s, m]
    flux_tau_down = conf.F_tau[s, m_prev]
    dE_temporal = 2.0 * params.tau_coupling * (flux_tau_up + flux_tau_down)
    return dE_spatial + dE_temporal
end

function flip_spatial!(conf::SimConfig, params::SimParams, s::Int, m::Int)
    @inbounds conf.ising_config_G[s, m] *= -1
    @inbounds for p_idx in conf.edge_to_flux[s]
        conf.F_G[p_idx, m] *= -1
    end
    m_prev = prev_time(m, params.M)
    @inbounds conf.F_tau[s, m] *= -1
    @inbounds conf.F_tau[s, m_prev] *= -1
end

function calc_dE_temporal(conf::SimConfig, params::SimParams, v::Int, m::Int)
    sum_flux_tau = 0.0    
    @inbounds for e_idx in conf.vert_to_edge[v]
        sum_flux_tau += conf.F_tau[e_idx, m]
    end    
    return 2.0 * params.tau_coupling * sum_flux_tau
end

function flip_temporal!(conf::SimConfig, params::SimParams, v::Int, m::Int)
    @inbounds conf.ising_config_tau[v, m] *= -1
    @inbounds for e_idx in conf.vert_to_edge[v]
        conf.F_tau[e_idx, m] *= -1
    end
end

function attempt_spin_flip!(conf::SimConfig, params::SimParams, s::Int, m::Int, G_or_tau::Int)    
    if G_or_tau == 1
        dE = calc_dE_spatial(conf, params, s, m)        
        if dE <= 0 || rand() < exp(-dE)
            flip_spatial!(conf, params, s, m)
            return true
        end
    else
        dE = calc_dE_temporal(conf, params, s, m)
        if dE <= 0 || rand() < exp(-dE)
            flip_temporal!(conf, params, s, m)
            return true
        end
    end
    return false
end

function sweep_move!(conf::SimConfig, params::SimParams)
    M = params.M
    nm_G = length(conf.edge_to_vert)
    nm_tau = length(conf.vert_to_edge)
    total_steps = (nm_G + nm_tau) * M
    prob_spatial = nm_G / (nm_G + nm_tau)
    for i=1:total_steps
        G_or_tau = (rand(params.rng) < prob_spatial) ? 1 : 0
        if G_or_tau == 1
            s = rand(params.rng, 1:nm_G)
            m = rand(params.rng, 1:M)
            attempt_spin_flip!(conf,params,s,m,G_or_tau)
        else
            s = rand(params.rng, 1:nm_tau)
            m = rand(params.rng, 1:M)
            attempt_spin_flip!(conf,params,s,m,G_or_tau)
        end
    end
end

function temporal_sweep_vectorized!(conf::SimConfig, params::SimParams)
    nm_tau = length(conf.vert_to_edge)
    M = params.M
    coupling_factor = 2.0 * params.tau_coupling
    dEs = zeros(Float64, M)
    randoms = zeros(Float64, M)
    
    for v in 1:nm_tau
        edges = conf.vert_to_edge[v]
        fill!(dEs, 0.0)
        for e_idx in edges
            @turbo for m in 1:M
                dEs[m] += conf.F_tau[e_idx, m]
            end
        end
        rand!(params.rng, randoms)
        @inbounds for m in 1:M
            delta_E = coupling_factor * dEs[m]
            if delta_E <= 0 || randoms[m] < exp(-delta_E)
                conf.ising_config_tau[v, m] *= -1
                for e_idx in edges
                    conf.F_tau[e_idx, m] *= -1
                end
            end
        end
    end
end

function spatial_sweep_vectorized!(conf::SimConfig, params::SimParams)
    nm_G = length(conf.edge_to_vert)
    M = params.M
    G_coupling_factor = 2.0 * params.G_coupling
    tau_coupling_factor = 2.0 * params.tau_coupling
    dEs = conf.G_dEs
    randoms = conf.G_randoms
    for s in 1:nm_G
        spatial_neighbors = conf.edge_to_flux[s]
        fill!(dEs, 0.0)        
        for p_idx in spatial_neighbors
            @turbo for m in 1:M
                dEs[m] += conf.F_G[p_idx, m]
            end
        end
        @turbo for m in 1:M
            dEs[m] *= G_coupling_factor
        end
        rand!(params.rng, randoms)  
              
        @inbounds for m in 1:2:M
            m_prev = (m == 1) ? M : m - 1
            flux_sum = conf.F_tau[s, m] + conf.F_tau[s, m_prev]            
            current_dE = dEs[m] + (tau_coupling_factor * flux_sum)
            if current_dE <= 0 || randoms[m] < exp(-current_dE)
                conf.ising_config_G[s, m] *= -1
                for p_idx in spatial_neighbors
                    conf.F_G[p_idx, m] *= -1
                end
                conf.F_tau[s, m] *= -1
                conf.F_tau[s, m_prev] *= -1
            end
        end
        
        @inbounds for m in 2:2:M
            m_prev = m - 1
            flux_sum = conf.F_tau[s, m] + conf.F_tau[s, m_prev]
            current_dE = dEs[m] + (tau_coupling_factor * flux_sum)
            if current_dE <= 0 || randoms[m] < exp(-current_dE)
                conf.ising_config_G[s, m] *= -1
                for p_idx in spatial_neighbors
                    conf.F_G[p_idx, m] *= -1
                end
                conf.F_tau[s, m] *= -1
                conf.F_tau[s, m_prev] *= -1
            end
        end
    end
end

function sweep_move_vectorized!(conf::SimConfig, params::SimParams)
    temporal_sweep_vectorized!(conf,params)
    spatial_sweep_vectorized!(conf,params)
end