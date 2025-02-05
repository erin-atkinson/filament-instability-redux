# Buoyancy flux
@inline function BFLUX_func(i, j, k, grid, w, b, mean_fields)
    b_dfm = mean_fields.b_dfm
    w_dfm = mean_fields.w_dfm
    return ℑzᵃᵃᶜ(i, j, k, grid, j_mean, w′b′ᶜᶜᶠ_func, w, w_dfm, b, b_dfm)
end

# Lateral and vertical shear production
@inline function LSP_func(i, j, k, grid, u, v, w, mean_fields)
    
    u_dfm = mean_fields.u_dfm
    v_dfm = mean_fields.v_dfm
    w_dfm = mean_fields.w_dfm
    
    x = ℑxᶜᵃᵃ(i, j, k, grid, j_mean, u′u′ᶠᶜᶜ_func, u, u_dfm) * ∂xᶜᶜᶜ(i, j, k, grid, u_dfm)
    y = ℑxᶜᵃᵃ(i, j, k, grid, j_mean, u′v′ᶠᶜᶜ_func, u, u_dfm, v, v_dfm) * ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, v_dfm)
    z = ℑxᶜᵃᵃ(i, j, k, grid, j_mean, u′w′ᶠᶜᶜ_func, u, u_dfm, w, w_dfm) * ℑxzᶜᵃᶜ(i, j, k, grid, ∂xᶠᶜᶠ, w_dfm)
    
    return -(x + y + z)
end

@inline function VSP_func(i, j, k, grid, u, v, w, mean_fields)
    
    u_dfm = mean_fields.u_dfm
    v_dfm = mean_fields.v_dfm
    w_dfm = mean_fields.w_dfm
    
    x = ℑzᵃᵃᶜ(i, j, k, grid, j_mean, w′u′ᶜᶜᶠ_func, w, w_dfm, u, u_dfm) * ℑxzᶜᵃᶜ(i, j, k, grid, ∂zᶠᶜᶠ, u_dfm)
    y = ℑzᵃᵃᶜ(i, j, k, grid, j_mean, w′v′ᶜᶜᶠ_func, w, w_dfm, v, v_dfm) * ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, v_dfm)
    z = ℑzᵃᵃᶜ(i, j, k, grid, j_mean, w′w′ᶜᶜᶠ_func, w, w_dfm) * ∂zᶜᶜᶜ(i, j, k, grid, w_dfm)
    
    return -(x + y + z)
end

@inline function DISP_func(i, j, k, grid, u, v, w, ν)
    # The dissipation
    # 9 terms
    u_x = ∂xᶜᶜᶜ(i, j, k, grid, u)
    u_y = ℑxyᶜᶜᵃ(i, j, k, grid, ∂yᶠᶠᶜ, u)
    u_z = ℑxzᶜᵃᶜ(i, j, k, grid, ∂zᶠᶜᶠ, u)
    
    v_x = ℑxyᶜᶜᵃ(i, j, k, grid, ∂xᶠᶠᶜ, v)
    v_y = ∂yᶜᶜᶜ(i, j, k, grid, v)
    v_z = ℑyzᵃᶜᶜ(i, j, k, grid, ∂zᶜᶠᶠ, v)
    
    w_x = ℑxzᶜᵃᶜ(i, j, k, grid, ∂xᶠᶜᶠ, w)
    w_y = ℑyzᵃᶜᶜ(i, j, k, grid, ∂yᶜᶠᶠ, w)
    w_z = ∂zᶜᶜᶜ(i, j, k, grid, w)
    
    return @inbounds ν[i, j, k] * (
        u_x^2 + u_y^2 + u_z^2
        + v_x^2 + v_y^2 + v_z^2
        + w_x^2 + w_y^2 + w_z^2
    )
end
