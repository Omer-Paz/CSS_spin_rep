
mutable struct Measurements
    i::Int64
    avg_spatial_flux::Vector{Float64}
    flux_corr::Vector{Float64}
end

function init_measurements(nm_meas::Int64)
    i=1
    avg_flux = zeros(Float64, nm_meas)
    flux_corr = zeros(Float64, nm_meas)
    return Measurements(i,avg_flux,flux_corr) 
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

function measure_all!(meas::Measurements,conf::SimConfig)
    measure_tot_flux!(meas,conf)
    measure_flux_corr!(meas,conf)
    meas.i+=1
end

