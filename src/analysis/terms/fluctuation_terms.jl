# Terms involving turbulent fluxes

@inline function u′u′ᶠᶜᶜ_func(i, j, k, grid, u, u_dfm)
    f′g′(i, j, k, grid, u, u_dfm, u, u_dfm)
end

@inline function w′u′ᶜᶜᶠ_func(i, j, k, grid, w, w_dfm, u, u_dfm)
    f′Gg′(i, j, k, grid, w, w_dfm, ℑxzᶜᵃᶠ, u, u_dfm)
end

@inline function u′v′ᶠᶜᶜ_func(i, j, k, grid, u, u_dfm, v, v_dfm)
    f′Gg′(i, j, k, grid, u, u_dfm, ℑxyᶠᶜᵃ, v, v_dfm)
end

@inline function w′v′ᶜᶜᶠ_func(i, j, k, grid, w, w_dfm, v, v_dfm)
    f′Gg′(i, j, k, grid, w, w_dfm, ℑyzᵃᶜᶠ, v, v_dfm)
end

@inline function u′w′ᶠᶜᶜ_func(i, j, k, grid, u, u_dfm, w, w_dfm)
    f′Gg′(i, j, k, grid, u, u_dfm, ℑxzᶠᵃᶜ, w, w_dfm)
end

@inline function w′w′ᶜᶜᶠ_func(i, j, k, grid, w, w_dfm)
    f′g′(i, j, k, grid, w, w_dfm, w, w_dfm)
end

@inline function u′b′ᶠᶜᶜ_func(i, j, k, grid, u, u_dfm, b, b_dfm)
    f′Gg′(i, j, k, grid, u, u_dfm, ℑxᶠᵃᵃ, b, b_dfm)
end

@inline function w′b′ᶜᶜᶠ_func(i, j, k, grid, w, w_dfm, b, b_dfm)
    f′Gg′(i, j, k, grid, w, w_dfm, ℑzᵃᵃᶠ, b, b_dfm)
end
