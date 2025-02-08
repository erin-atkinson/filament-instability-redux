# hovmoller.jl
# Creates hovmoller plots of the mean buoyancy gradient
# and dissipation

using Statistics: mean

function hovmoller_figure(
        foldername;
        fig_kw=(; ),
        ax_b_kw=(; ),
        ax_ϵ_kw=(; ),
        ht_b_kw=(; ),
        ht_ϵ_kw=(; ),
        ct_b_kw=(; ),
        z_top=0,
        z_bottom=-0.01,
        σ=0,
        marked_times=[]
    )
    
    L, T = scales()
    
    fig_kw = (; 
        size=(750, 230), 
        fig_kw...
    )
    fig = Figure(; fig_kw...)
    
    sp = simulation_parameters(foldername)
    iterations, ts = iterations_times(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    
    ax_b_kw = (;
        xlabel=L"ft / 2\pi",
        ylabel=L"x / \text{km}",
        limits=(0, 4, -2sp.L * L / 1000, 2sp.L * L / 1000),
        ax_b_kw...
    )
    
    ax_ϵ_kw = (;
        xlabel=L"ft / 2\pi",
        ylabel=L"x / \text{km}",
        limits=(0, 4, -2sp.L * L / 1000, 2sp.L * L / 1000),
        ax_ϵ_kw...
    )
    
    ax_b = Axis(fig[1, 1]; ax_b_kw...)
    
    ax_ϵ = Axis(fig[1, 3]; ax_ϵ_kw...)
    hideydecorations!(ax_ϵ; ticks=false)
    
    bottom_index, top_index = zᶜbounds(foldername, z_bottom, z_top)
    reduce_func(field) = mean(field[:, bottom_index:top_index]; dims=2)[:, 1]
    
    DFM = joinpath(foldername, "DFM.jld2")
    DISP = joinpath(foldername, "DISP.jld2")
    
    b = filt(timeseries_of(reduce_func, DFM, "b_dfm", iterations), σ)
    ϵ = filt(timeseries_of(reduce_func, DISP, "DISP_dfm", iterations) .* (1e7L^2 / T^3), σ)
    b_x = filt(timeseries_of(reduce_func, DISP, "b_x", iterations) .* (1e7 / T^2), σ)
    b_x_colors = to_colormap(:balance)
    ϵ_colors = reverse(to_colormap(:acton100))
    
    ht_b_kw = (;
        colormap=b_x_colors,
        colorrange=(-4, 4),
        lowclip=b_x_colors[1],
        highclip=b_x_colors[end],
        ht_b_kw...
    )
    
    ht_ϵ_kw = (;
        colormap=ϵ_colors,
        colorrange=(0, 1),
        highclip=ϵ_colors[end],
        ht_ϵ_kw...
    )
    
    ht_b = heatmap!(ax_b, sp.f*ts/2π, xsᶜ*L/1000, b_x; ht_b_kw...)
    ht_ϵ = heatmap!(ax_ϵ, sp.f*ts/2π, xsᶜ*L/1000, ϵ; ht_ϵ_kw...)

    ct_b_kw = (;
        levels=range(-300, -227, 12),
        color=(:black, 0.3),
        linewidth=2,
        ct_b_kw...
    )
    
    contour!(ax_b, sp.f*ts/2π, xsᶜ*L/1000, b; ct_b_kw...)
    contour!(ax_ϵ, sp.f*ts/2π, xsᶜ*L/1000, b; ct_b_kw...)
    
    
    Colorbar(fig[1, 2], ht_b, label=L"10^{7}\overline{b}_{,x} / \text{s}^{-2}", ticks=WilkinsonTicks(5; k_max=5))
    Colorbar(fig[1, 4], ht_ϵ, label=L"10^{7}\overline{\varepsilon} / \text{m}^2\text{s}^{-3}", ticks=WilkinsonTicks(5; k_max=6))
    
    subfig_label!(fig[1, 1], 1)
    subfig_label!(fig[1, 3], 2)
    
    for t in marked_times
        x = [t, t]
        y = [ax_b_kw.limits[3], ax_b_kw.limits[4]]
        lines!(ax_b, x, y; color=(:red, 0.5), linestyle=:dash)
        y = [ax_ϵ_kw.limits[3], ax_ϵ_kw.limits[4]]
        lines!(ax_ϵ, x, y; color=(:red, 0.5), linestyle=:dash)
    end
    
    return fig
end