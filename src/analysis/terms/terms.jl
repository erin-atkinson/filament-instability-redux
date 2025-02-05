using Oceananigans
using Oceananigans.Operators

# Anytime I need to calculate something that I will probably reuse it goes here

@inline dfm(a) = Field(Average(a; dims=2))

@inline function j_mean(i, j, k, grid, f, args...)
    total = mapreduce(+, (grid.Hy+1):(grid.Hy+grid.Ny)) do _j
        f(i, _j, k, grid, args...)
    end
    total / grid.Ny
end

@inline fg(i, j, k, grid, f, g) = @inbounds f[i, j, k] * g[i, j, k]
@inline fGg(i, j, k, grid, f, G, g) = @inbounds f[i, j, k] * G(i, j, k, grid, g)
@inline fG1G2g(i, j, k, grid, f, G1, G2, g) = @inbounds f[i, j, k] * G1(i, j, k, grid, G2, g)

@inline f′(i, j, k, grid, f, f_dfm) = @inbounds f[i, j, k] - f_dfm[i, 1, k]
@inline f′g′(i, j, k, grid, f, f_dfm, g, g_dfm) = f′(i, j, k, grid, f, f_dfm) * f′(i, j, k, grid, g, g_dfm)
@inline f′Gg′(i, j, k, grid, f, f_dfm, G, g, g_dfm) = f′(i, j, k, grid, f, f_dfm) * G(i, j, k, grid, f′, g, g_dfm)


# Advection by the mean and deformation flow
@inline function u_adv_func(i, j, k, grid, mean_fields, args...)
    return -ℑxᶜᵃᵃ(i, j, k, grid, fGg, mean_fields.u_dfm, ∂xᶠᶜᶜ, args...)
end

@inline function w_adv_func(i, j, k, grid, mean_fields, args...)
    return -ℑzᵃᵃᶜ(i, j, k, grid, fGg, mean_fields.w_dfm, ∂zᶜᶜᶠ, args...)
end
