#=
parameters.jl
    Input parameters are non-dimensional numbers, here any dependent parameters are derived
    Ro=1: Maximum magnitude of Rossby number
    Ri=0.6: Minimum Richardson number
    Ek=nothing: (Turbulent) Ekman number, for closure. If nothing then a Smagorinsky-Lilly default is used
    Pr=1: Prandtl number for buoyancy diffusion
    α=0.25: Front width - separation ratio
    λ=0.05: Fractional width of transition to deep water
    δ=-0.25: Fractional height change of transition to deep water across filament

    λ ≪ δ ≪ 1 is required for the correct computation of parameters. (factors of two is enough)
=#
using SpecialFunctions
# Filament shape function
γ(s, δ, α) = -1 + (δ / 2) * (erf(s + 1/(2α)) - erf(s - 1/(2α)))
# Function to maximise for Ro
square_curvature(x, δ, α) = -2*((γ(x, δ, α)^2 - (γ(x+1e-5, δ, α)^2 + γ(x-1e-5, δ, α)^2) / 2) / 1e-10)
# Function to maximise for Ri
grad_squared(x, δ, α) = ((γ(x+1e-5, δ, α) - γ(x-1e-5, δ, α)) / (2e-5))^2

default_inputs = (; 
    Ro=1, Ri=0.6, α=0.25, β=0.1, Nx=128, Ny=128, Nz=128,
    Q=0.0, noise_amplitude=0.0, init_time=0.0, init_rate=0.0,
    preinit=false, preinit_path=""
)

# Q=5e-5 corresponds to about 100Wm-2

@inline function create_simulation_parameters(input_parameters=(; ))
    ip = (; default_inputs..., input_parameters...)
    let Ro=ip.Ro, Ri=ip.Ri, α=ip.α, β=ip.β, Nx=ip.Nx, Ny=ip.Ny, Nz=ip.Nz,
        Q=ip.Q, noise_amplitude=ip.noise_amplitude, init_time=ip.init_time, init_rate=ip.init_rate,
        preinit=ip.preinit, preinit_path = ip.preinit_path
        # Setting variables
        # The length and time scales are set as the distance between fronts
        # and the coriolis frequency
        # Distance between fronts
        L = 1
        # Coriolis frequency
        f = 1
        # Thermocline transition width
        λ = 0.05
        # MLD change
        δ = -0.25
        
        # Derived variables
        # MLD
        H = β*L
        δH = δ * H
        
        # Front width
        ℓ = α * L
        
        # Difference in stratification from Rossby number
        ΔN² = -(2Ro * f^2 * ℓ^2) / (H^2 * abs(maximum([square_curvature(x, δ, α) for x in range(-3L, 3L, 1000)])))
        # Boundary layer stratification to reach a specific min Ri
        Nb² = ((Ri * ΔN²^2 * H^2) / (f^2 * ℓ^2)) * maximum([grad_squared(x, δ, α) for x in range(-3L, 3L, 1000)])
        
        # Convenience
        Nb = sqrt(Nb²)
        N₀² = Nb² - ΔN²
        N₀ = sqrt(N₀²)
        
        # Domain size
        Lx = 10L
        Ly = (Ny / Nx) * Lx
        Lz = 2.5H
        
        # Damping constant
        c = 0.5
        
        # Remove the filament...
        preinit && (δ = 0)
        (; 
            Ro, Ri, α, β, L, H, f, λ, δ, δH, ℓ, ΔN², Nb², N₀²,
            Nb, N₀, Lx, Ly, Lz, Nx, Ny, Nz, c, Q, noise_amplitude,
            init_time, init_rate, 
            preinit, preinit_path
        )
    end
end

@inline function create_simulation_parameters(; input_parameters...)
    create_simulation_parameters(input_parameters)
end
    