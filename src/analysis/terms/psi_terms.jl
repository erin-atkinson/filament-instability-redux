
# Linear terms
function Lψ_func(i, j, k, grid, v₀, b₀, ψ, sp)
    
    ζ₀ = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, v₀)
    S₀ = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, v₀)
    M₀² = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b₀)
    N₀² = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, b₀)
    
    ψ_zz = ∂zᶜᶜᶜ(i, j, k, grid, ∂zᶜᶜᶠ, ψ)
    ψ_xz = ℑxzᶜᵃᶜ(i, j, k, grid, ∂xᶠᶜᶠ, ∂zᶜᶜᶠ, ψ)
    ψ_xx = ∂xᶜᶜᶜ(i, j, k, grid, ∂xᶠᶜᶜ, ψ)
    
    return @inbounds sp.f^2 * ψ_zz[i, j, k] + sp.f * ζ₀[i, j, k] * ψ_zz[i, j, k] - sp.f * S₀[i, j, k] * ψ_xz[i, j, k] - M₀²[i, j, k] * ψ_xz[i, j, k] + N₀²[i, j, k] * ψ_xx[i, j, k]
end
function ∇²ψ_func(i, j, k, grid, ψ)
    return ∂xᶜᶜᶜ(i, j, k, grid, ∂xᶠᶜᶜ, ψ) + ∂zᶜᶜᶜ(i, j, k, grid, ∂zᶜᶜᶠ, ψ)
end

# Terms induced by turbulent fluxes
@inline function Gh_t_func(i, j, k, grid, fields, mean_fields, sp)
    
    u = fields.u
    #v = fields.v
    w = fields.w
    #b = fields.b
    
    u_dfm = mean_fields.u_dfm
    #v_dfm = mean_fields.v_dfm
    w_dfm = mean_fields.w_dfm
    #b_dfm = mean_fields.b_dfm
    
    A = -ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, ∂xᶜᶜᶜ, j_mean, u′u′ᶠᶜᶜ_func, u, u_dfm)
    B = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂xᶜᶜᶜ, j_mean, u′w′ᶠᶜᶜ_func, u, u_dfm, w, w_dfm)
    
    return A + B
end

@inline function Gh_nt_func(i, j, k, grid, fields, mean_fields, sp)
    
    u = fields.u
    v = fields.v
    #w = fields.w
    b = fields.b
    
    u_dfm = mean_fields.u_dfm
    v_dfm = mean_fields.v_dfm
    #w_dfm = mean_fields.w_dfm
    b_dfm = mean_fields.b_dfm
    
    A = -sp.f * ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, ∂xᶜᶜᶜ, j_mean, u′v′ᶠᶜᶜ_func, u, u_dfm, v, v_dfm)
    B = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂xᶜᶜᶜ, j_mean, u′b′ᶠᶜᶜ_func, u, u_dfm, b, b_dfm)
    
    return A + B
end

@inline function Gv_t_func(i, j, k, grid, fields, mean_fields, sp)
    
    u = fields.u
    #v = fields.v
    w = fields.w
    #b = fields.b
    
    u_dfm = mean_fields.u_dfm
    #v_dfm = mean_fields.v_dfm
    w_dfm = mean_fields.w_dfm
    #b_dfm = mean_fields.b_dfm
    
    A = -ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, ∂zᶜᶜᶜ, j_mean, w′u′ᶜᶜᶠ_func, w, w_dfm, u, u_dfm)
    B = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂zᶜᶜᶜ, j_mean, w′w′ᶜᶜᶠ_func, w, w_dfm)
    
    return A + B
end

@inline function Gv_nt_func(i, j, k, grid, fields, mean_fields, sp)
    
    #u = fields.u 
    v = fields.v
    w = fields.w
    b = fields.b
    
    #u_dfm = mean_fields.u_dfm
    v_dfm = mean_fields.v_dfm
    w_dfm = mean_fields.w_dfm
    b_dfm = mean_fields.b_dfm
    
    A = -sp.f * ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, ∂zᶜᶜᶜ, j_mean, w′v′ᶜᶜᶠ_func, w, w_dfm, v, v_dfm)
    B = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂zᶜᶜᶜ, j_mean, w′b′ᶜᶜᶠ_func, w, w_dfm, b, b_dfm)
    
    return A + B
end