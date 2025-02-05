#

function w_slices_figure(
        foldername,
        t1, 
        t2,
        h; 
        fig_kw=(; ), 
        ax_w_kw=(; ), 
        ax_wh_kw=(; ),
        ht_w_kw=(; ),
        ht_wh_kw=(; ),
        ct_b_kw=(; ),
        σ=0
    )
    
    L, T = scales()
    b_levels = range(-300, -227, 16)
    w_colors = to_colormap(:balance)
    
    fig_kw = (;
        size=(750, 500),
        fig_kw...
    )
    fig = Figure(; fig_kw...)
    
    sp = simulation_parameters(foldername)
    iterations, ts = iterations_times(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    
    z_ind = zᶠbounds(foldername, h)
    
    i1 = iterations[tbounds(foldername, t1)]
    i2 = iterations[tbounds(foldername, t2)]
    
    OUTPUT = joinpath(foldername, "output.jld2")
    
    w1 = filt(get_field(a->a[:, 1, :], OUTPUT, "w", i1) * 100L / T, σ, 0)
    w2 = filt(get_field(a->a[:, 1, :], OUTPUT, "w", i2) * 100L / T, σ, 0)
    wh = transpose(filt(get_field(a->a[:, :, z_ind], OUTPUT, "w", i2) * 100L / T, σ))
    
    b1 = filt(get_field(a->a[:, 1, :], OUTPUT, "b", i1), σ, 0)
    b2 = filt(get_field(a->a[:, 1, :], OUTPUT, "b", i2), σ, 0)
    bh = transpose(filt(get_field(a->a[:, :, z_ind], OUTPUT, "b", i2), σ))
    
    ax_w_kw = (;
        xlabel=L"x / \text{km}",
        ylabel=L"z / \text{m}", 
        limits=(-2.2sp.L * L/1000, -0.5sp.L * L/1000, -1.2sp.H * L, 0*L),
        xlabelsize=16,
        ylabelsize=16,
        yticks=WilkinsonTicks(3; k_max=5),
        ax_w_kw...
    )
    
    ax_wh_kw = (;
        xlabel=L"y / \text{km}", 
        ylabel=L"x / \text{km}", 
        limits=(-5sp.L * L/1000, 5sp.L * L/1000, -2.2sp.L * L/1000, -0.5sp.L * L/1000), 
        xlabelsize=16, 
        ylabelsize=16,
        ax_wh_kw...
    )
    
    ax_w1 = Axis(fig[1, 1]; ax_w_kw...)
    ax_w2 = Axis(fig[1, 2]; ax_w_kw...)
    ax_wh = Axis(fig[2, 1:2]; ax_wh_kw...)
    
    hideydecorations!.(ax_w2; ticks=false)
    
    w_max = maximum(abs, w2)
    
    ht_w_kw = (;
        colormap=w_colors,
        colorrange=(-w_max, w_max),
        lowclip=w_colors[1],
        highclip=w_colors[end],
        ht_w_kw...
    )
    
    ht_wh_kw = (;
        colormap=w_colors,
        colorrange=(-w_max, w_max),
        lowclip=w_colors[1],
        highclip=w_colors[end],
        ht_wh_kw...
    )
    
    ct_b_kw=(; 
        color=(:black, 0.5),
        levels=b_levels,
        linewidth=2,
        ct_b_kw...
    )
    
    ht_w1 = heatmap!(ax_w1, xsᶜ*L/1000, zsᶠ*L, w1; ht_w_kw...)
    ht_w2 = heatmap!(ax_w2, xsᶜ*L/1000, zsᶠ*L, w2; ht_w_kw...)
    ht_wh = heatmap!(ax_wh, ysᶜ*L/1000, xsᶜ*L/1000, wh; ht_wh_kw...)
    
    contour!(ax_w1, xsᶜ*L/1000, zsᶜ*L, b1; ct_b_kw...)
    contour!(ax_w2, xsᶜ*L/1000, zsᶜ*L, b2; ct_b_kw...)
    contour!(ax_wh, ysᶜ*L/1000, xsᶜ*L/1000, bh; ct_b_kw...)
    
    lines!(ax_w2, [ax_w_kw.limits[1], ax_w_kw.limits[2]], [h, h]*L; color=(:red, 0.5), linestyle=:dash)
    
    Colorbar(fig[1:2, 3], ht_w2; label=L"w / \text{cms}^{-1}")
    
    colgap!(fig.layout, 1, 30)
    #colgap!(fig.layout, 2, 30)
    
    subfig_label!(fig[1, 1], 1)
    subfig_label!(fig[1, 2], 2)
    subfig_label!(fig[2, 1:2], 3)
    
    fig
end