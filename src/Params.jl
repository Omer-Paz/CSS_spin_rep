using Random

struct SimParams
    # == Physical Parameters ==
    beta::Float64          # Inverse temperature
    M::Int                 # Number of Trotter slices
    Jz::Float64            # Spatial term coefficient 
    h::Float64             # Transverse field coefficient
    epsilon::Float64       # Perturbation parameter / or dtau context
    tau_coupling::Float64  # Coupling in time direction
    G_coupling::Float64    # Coupling in spatial direction
    
    # == Monte Carlo Parameters ==
    nm_sweep::Int          # Sweeps between measurements
    nm_therm::Int          # Sweeps for thermalization
    nm_meas::Int           # Number of measurements
    
    # == RNG ==
    rng::MersenneTwister

    # Constructor
    function SimParams(beta, M, Jz, h, epsilon, nm_sweep, nm_therm, nm_meas; seed=nothing)
        dtau = beta / M
        curr_G_coupling = dtau * Jz 
        curr_tau_coupling = -0.5 * log(tanh(dtau * h))        
        final_seed = isnothing(seed) ? rand(RandomDevice(), UInt64) : seed
        rng = MersenneTwister(final_seed)
        new(beta, M, Jz, h, epsilon, 
            curr_tau_coupling, curr_G_coupling, 
            nm_sweep, nm_therm, nm_meas, rng)
    end
end