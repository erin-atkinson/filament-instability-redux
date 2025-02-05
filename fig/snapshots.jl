# snapshots.jl
# Snapshots of the simulation at some times
# Down-front velocity and streamfunction

function snapshots_figure(
        foldername, 
        times;
        fig_kw=(; ),
        ax_kw=(; ),
        ht_kw=(; ),
        ct_b_kw=(; ),
        ct_ψ_kw=(; ),
        σ=0,
        nψ=5,
    )
    
    L, T = scales()
    
    n_plots = length(times)
    
    fig_kw = (; 
        size=(750, 230), 
        fig_kw...
    )
    fig = Figure(; fig_kw...)
    
    sp = simulation_parameters(foldername)
    iterations, ts = iterations_times(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    
    iterations = map(t->iterations[argmin(abs.(t .- ts))], times)
    
    ax_kw = (;
        xlabel=L"x / \text{km}",
        ylabel=L"z / \text{m}",
        yticks=WilkinsonTicks(5; k_max=5),
        title=L"",
        limits=(-3.3sp.L*L/1000, 3.3sp.L*L/1000, -1.2sp.H*L, 0), 
        xlabelsize=16, 
        ylabelsize=16,
        ax_kw...
    )
    
    # Get the data
    DFM = joinpath(foldername, "DFM.jld2")
    vs = map(iterations) do iteration
        filt(get_field(DFM, "v_dfm", iteration), σ) * L / T
    end
    bs = map(iterations) do iteration
        filt(get_field(DFM, "b_dfm", iteration), σ)
    end
    ψs = map(iterations) do iteration
        filt(get_field(DFM, "ψ", iteration), σ)
    end
    
    # Aim for three contours each sign
    ψ_max = mapreduce(max, ψs) do ψ
        maximum(abs, ψ)
    end
    ψ_levels = range(ψ_max / (nψ + 0.5), ψ_max, nψ)
    
    # Make all the axes
    axes = map(enumerate(times)) do (i, t)
        title = join([ax_kw.title, L"ft/2\pi=%$(round(sp.f*t / 2π; digits=1))"], L"\quad")
        title = L"%$title"
        ax = Axis(fig[1, i]; ax_kw..., title)
        i > 1 && hideydecorations!(ax; ticks=false)
        ax
    end
    
    # Plot v heatmaps
    v_colors = to_colormap(:balance)
    ht_kw = (;
        colormap=v_colors,
        colorrange=(-0.2, 0.2),
        lowclip=v_colors[1],
        highclip=v_colors[end],
        ht_kw...
    )
    ht_vs = map(axes, vs) do ax, v
        heatmap!(ax, xsᶜ*L/1000, zsᶜ*L, v; ht_kw...)
    end
    
    # Plot b contours
    ct_b_kw = (;
        levels=range(-300, -227, 12),
        color=(:black, 0.5),
        linewidth=2,
        ct_b_kw...
    )
    ct_bs = map(axes, bs) do ax, b
        contour!(ax, xsᶜ*L/1000, zsᶜ*L, b; ct_b_kw...)
    end
    
    # Plot streamfunction
    ct_ψ_kw = (;
        linewidth=1.5,
        linestyle=:dash,
        colormap=reverse(to_colormap(:PuOr)),
        ct_ψ_kw...
    )
    ct_ψs = map(axes, ψs) do ax, ψ
        ct1 = contour!(ax, xsᶜ*L/1000, zsᶜ*L, ψ; ct_ψ_kw..., levels=-ψ_levels, color=ct_ψ_kw.colormap[1])
        ct2 = contour!(ax, xsᶜ*L/1000, zsᶜ*L, ψ; ct_ψ_kw..., levels=ψ_levels, color=ct_ψ_kw.colormap[end])
        [ct1, ct2]
    end
    
    Colorbar(fig[1, n_plots+1], ht_vs[1], label=L"\overline{v} / \text{ms}^{-1}", ticks=WilkinsonTicks(5; k_max=5))
    
    for i in 1:n_plots
        subfig_label!(fig[1, i], i)
    end
    
    return fig
end