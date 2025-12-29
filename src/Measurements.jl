using LoopVectorization

mutable struct Measurements
    i::Int64
    avg_spatial_flux::Vector{Float64}
    flux_corr::Vector{Float64}
    wilson_loop::Matrix{Float64}
end

function init_measurements(nm_meas::Int64,N_spatial::Int64)
    i=1
    avg_flux = zeros(Float64, nm_meas)
    flux_corr = zeros(Float64, nm_meas)
    wilson_loop = zeros(Float64, N_spatial, nm_meas)
    return Measurements(i,avg_flux,flux_corr,wilson_loop) 
end

function measure_flux_corr!(meas::Measurements,conf::SimConfig)
    M = size(conf.F_G,2)
    for m=1:M
        meas.flux_corr[meas.i] += conf.F_G[1,m]*conf.F_G[6,m]
    end
    meas.flux_corr[meas.i] /= M
end

function measure_tot_flux!(meas::Measurements,conf::SimConfig)
    meas.avg_spatial_flux[meas.i] = sum(Float64, conf.F_G) / length(conf.F_G)
end

function measure_wilson_loop!(meas::Measurements,conf::SimConfig)
    N_verts, M = size(conf.ising_config_tau)
    polyakov_lines = ones(Int8, N_verts)
    @turbo for m in 1:M
        for v in 1:N_verts
            polyakov_lines[v] *= conf.ising_config_tau[v, m]
        end
    end
    p1 = polyakov_lines[1]
    idx = meas.i
    @turbo for s in 1:N_verts
        meas.wilson_loop[s, idx] = p1 * polyakov_lines[s]
    end
end

function measure_all!(meas::Measurements,conf::SimConfig)
    measure_tot_flux!(meas,conf)
    measure_flux_corr!(meas,conf)
    measure_wilson_loop!(meas,conf)
    meas.i+=1
end

