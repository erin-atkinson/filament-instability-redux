@inline function η_x_func(i, j, k, grid, velocities)
    return ∂yᶜᶜᶜ(i, j, k, grid, ℑyzᵃᶠᶜ, velocities.w) - ∂zᶜᶜᶜ(i, j, k, grid, ℑyzᵃᶜᶠ, velocities.v)
end

@inline function η_y_func(i, j, k, grid, velocities)
    return ∂zᶜᶜᶜ(i, j, k, grid, ℑxzᶜᵃᶠ, velocities.u) - ∂xᶜᶜᶜ(i, j, k, grid, ℑxzᶠᵃᶜ, velocities.w)
end

@inline function η_z_func(i, j, k, grid, velocities)
    return ∂xᶜᶜᶜ(i, j, k, grid, ℑxyᶠᶜᵃ, velocities.v) - ∂yᶜᶜᶜ(i, j, k, grid, ℑxyᶜᶠᵃ, velocities.u)
end

@inline function full_q_func(i, j, k, grid, u, v, w, b, f)
    # Vorticity
    a = (∂yᶜᶜᶜ(i, j, k, grid, ℑyzᵃᶠᶜ, w) - ∂zᶜᶜᶜ(i, j, k, grid, ℑyzᵃᶜᶠ, v)) * ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b)
    b = (∂zᶜᶜᶜ(i, j, k, grid, ℑxzᶜᵃᶠ, u) - ∂xᶜᶜᶜ(i, j, k, grid, ℑxzᶠᵃᶜ, w)) * ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, b)
    c = (∂xᶜᶜᶜ(i, j, k, grid, ℑxyᶠᶜᵃ, v) - ∂yᶜᶜᶜ(i, j, k, grid, ℑxyᶜᶠᵃ, u) + f) * ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, b)
    
    return a + b + c
end

@inline function q_func(i, j, k, grid, v_dfm, b_dfm, f)
    
    a = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, v_dfm)
    b = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b_dfm)
    
    c = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, v_dfm)
    d = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, b_dfm)
    
    return d * (a + f) - b * c
end

@inline function UADV_func(i, j, k, grid, u, q)
    return -ℑxᶜᵃᵃ(i, j, k, grid, fGg, u, ∂xᶠᶜᶜ, q)
end

@inline function WADV_func(i, j, k, grid, w, q)
    return -ℑzᵃᵃᶜ(i, j, k, grid, fGg, w, ∂zᶜᶜᶠ, q)
end
