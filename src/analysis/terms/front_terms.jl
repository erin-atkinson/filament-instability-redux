@inline function T_adv_func(i, j, k, grid, mean_fields)
    
    u_dfm = mean_fields.u_dfm
    w_dfm = mean_fields.w_dfm
    b_dfm = mean_fields.b_dfm
    
    u_x = ∂xᶜᶜᶜ(i, j, k, grid, u_dfm)
    w_x = ℑxzᶜᵃᶜ(i, j, k, grid, ∂xᶠᶜᶠ, w_dfm)
    
    b_x = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b_dfm)
    b_z = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, b_dfm)
    
    return -(u_x * b_x + w_x * b_z) * ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b_dfm)
end

@inline function ∇b²_func(i, j, k, grid, mean_fields)
    return 0.5 * ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, mean_fields.b_dfm)^2
end

@inline function ∇b²_adv_func(i, j, k, grid, mean_fields, ∇b²)
    return -(u_adv_func(i, j, k, grid, mean_fields, ∇b²) + w_adv_func(i, j, k, grid, mean_fields, ∇b²))
end


@inline function Fv_func(i, j, k, grid, fields, mean_fields)
    
    w = fields.w
    b = fields.b
    
    w_dfm = mean_fields.w_dfm
    b_dfm = mean_fields.b_dfm
    
    return -ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂zᶜᶜᶜ, j_mean, w′b′ᶜᶜᶠ_func, w, w_dfm, b, b_dfm) * ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b_dfm)
end

@inline function Fh_func(i, j, k, grid, fields, mean_fields)
    
    u = fields.u
    b = fields.b
    
    u_dfm = mean_fields.u_dfm
    b_dfm = mean_fields.b_dfm
    
    return -ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂xᶜᶜᶜ, j_mean, u′b′ᶠᶜᶜ_func, u, u_dfm, b, b_dfm) * ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b_dfm)
end