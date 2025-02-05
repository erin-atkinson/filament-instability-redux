#=
forcing.jl
    Create a forcing function that relaxes the fields to zero in the bottom H of the domain
    at a rate c * internal wave freq that decays quadratically away from the boundary
=#
using Oceananigans.Operators

# Damping mask profile
σ₁(x) = abs(x) > 1 ? 0 : (1 - abs(x))^2
mask(x, y, z) = σ₁((z + sp.Lz) / sp.H)

rate = sp.c * sp.N₀ / 2π

# Optional damping to reference b if this is specified?

@inline function init_forcing(i, j, k, grid, clock, model_fields, p)
    clock.time >= 0 && return 0
    return @inbounds -sp.init_rate * (model_fields.b[i, j, k] - p.b₀[i, j, k])
end

b_op = KernelFunctionOperation{Center, Nothing, Center}(grid) do i, j, k, grid
    @inbounds b₀(grid.xᶜᵃᵃ[i], 0, grid.zᵃᵃᶜ[k])
end

b₀_field = Field(b_op)
compute!(b₀_field)

b_forcing = (
    Forcing(init_forcing; parameters=(; b₀=b₀_field), discrete_form=true),
    Relaxation(; rate, mask, target=(x, y, z, t)->b₀(x, y, z))
)
# Damping to the reference filament state at the bottom of the domain
forcing = (; 
    u=Relaxation(; rate, mask, target=(x, y, z, t)->0),
    v=Relaxation(; rate, mask, target=(x, y, z, t)->0),
    w=Relaxation(; rate, mask, target=(x, y, z, t)->0),
    b=b_forcing
)