#=
filament_state.jl
    Functions that define the nominal buoyancy state and associated thermal wind
=#
using SpecialFunctions: erf

# Transition function between the deep water and mixed layer
g(s) = (s + log(2*cosh(s)))/2
g′(s) = (1 + tanh(s)) / 2

# MLD function
γ(s, δ, α) = -1 + (δ / 2) * (erf(s + 1/(2α)) - erf(s - 1/(2α)))
∂γ∂s(s, δ, α) = (δ / sqrt(π)) * (exp(-(s + 1/(2α))^2) - exp(-(s - 1/(2α))^2))

# Scaled functions
filament_h(x) = sp.H*γ(x/sp.ℓ, sp.δ, sp.α)
filament_∂xh(x) = sp.H*∂γ∂s(x/sp.ℓ, sp.δ, sp.α)/sp.ℓ

# Velocity and buoyancy
b₀(x, y, z) = sp.N₀^2 * z + sp.H * sp.λ * (sp.Nb^2 - sp.N₀^2) * g((z - filament_h(x)) / (sp.λ * sp.H))
v₀(x, y, z) =  -(sp.H*sp.λ / sp.f) * filament_∂xh(x) * (sp.Nb^2 - sp.N₀^2) * g((z - filament_h(x)) / (sp.λ * sp.H))

# Noise in u and w?
u₀(x, y, z) = 2 * sp.noise_amplitude * (rand() - 0.5) * g((z + (1-4sp.λ) * sp.H) / (sp.λ * sp.H))
w₀(x, y, z) = 2 * sp.noise_amplitude * (rand() - 0.5) * g((z + (1-4sp.λ) * sp.H) / (sp.λ * sp.H))
