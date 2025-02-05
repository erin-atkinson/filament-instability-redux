# noise.jl
# Creates a figure with a slice of the initial vertical velocity
# and the horizontally-averaged turbulent kinetic energy

using Statistics: mean

@inline reduce_func(a) = mean(a.^2; dims=(1, 2))[1, 1, :] ./ 2

function noise_figure(
        foldername; 
        fig_kw=(; ), 
        ax_w_kw=(; ), 
        ax_tke_kw=(; ), 
        ht_kw=(; ), 
        ln_kw=(; )
    )
    
    L, T = scales()
    
    fig_kw = (; 
        size=(750, 230), 
        fig_kw...
    )
    
    sp = simulation_parameters(foldername)
    iterations, ts = iterations_times(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    
    ax_w_kw = (;
        xlabel=L"x / \text{km}",
        ylabel=L"z / \text{m}",
        limits=(-sp.L * L / 1000, sp.L * L / 1000, -sp.H * L, 0),
        ax_w_kw...
    )
    ax_tke_kw = (; 
        xlabel=L"10^{6}\text{KE}_\text{avg} / \text{m}^2\text{s}^{-2}",
        ylabel=L"z / \text{m}",
        limits=(nothing, nothing, -sp.H * L, 0),
        ax_tke_kw...
    )
    
    output_filename = joinpath(foldername, "output.jld2")
    
    tke = let 
        u = get_field(reduce_func, output_filename, "u", iterations[end])
        v = get_field(reduce_func, output_filename, "v", iterations[end])
        w = get_field(reduce_func, output_filename, "w", iterations[end])
        w = (w[1:end-1] .+ w[2:end]) ./ 2
        
        (u .+ v .+ w) .* (1e6 * L^2 / T^2)
    end
    
    w = get_field(w->w[:, 1, :] * (1000L / T), output_filename, "w", iterations[end])
    
    fig = Figure(; fig_kw...)
    
    ax_w = Axis(fig[1, 1:2]; ax_w_kw...)
    ax_tke = Axis(fig[1, 4]; ax_tke_kw...)
        
    # Left plot is a slice of the vertical velocity
    ht_kw = (;
        colormap=:balance,
        colorrange=(-maximum(abs, w), maximum(abs, w)),
        ht_kw...
    )
    ht_w = heatmap!(ax_w, xsᶜ * L / 1000, zsᶠ * L, w; ht_kw...)
    
    # Right plot is kinetic energy as a function of depth
    ln_kw = (;
        ln_kw...
    )
    ln_tke = lines!(ax_tke, tke, zsᶜ * L; ln_kw...)
    
    Colorbar(fig[1, 3], ht_w, label=L"w_\text{init} / \text{mm s}^{-1}")
    
    subfig_label!(fig[1, 1:2], 1)
    subfig_label!(fig[1, 4], 2)
    
    hideydecorations!(ax_tke; ticks=false, grid=false)
    
    fig
end