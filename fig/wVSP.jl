function wVSP_figure(
        foldername,
        t,
        h; 
        fig_kw=(; ), 
        ax_w_kw=(; ), 
        ax_vsp_kw=(; ), 
        ax_wh_kw=(; ),
        ht_w_kw=(; ),
        ht_vsp_kw=(; ),
        ht_wh_kw=(; ),
        ct_b_kw=(; ),
        σ=0,
        σh=0
    )
    
    L, T = scales()
    b_levels = range(-300, -227, 16)
    w_colors = to_colormap(:balance)
    vsp_colors = reverse(to_colormap(:curl))
    
    fig_kw = (;
        size=(750, 500),
        fig_kw...
    )
    
    fig = Figure(; fig_kw...)
    
    sp = simulation_parameters(foldername)
    iterations, ts = iterations_times(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    
    z_ind = zᶠbounds(foldername, h)
    
    i = iterations[tbounds(foldername, t)]
    
    OUTPUT = joinpath(foldername, "output.jld2")
    TKE = joinpath(foldername, "TKE.jld2")
    DFM = joinpath(foldername, "DFM.jld2")
    
    w = filt(get_field(a->a[:, 1, :], OUTPUT, "w", i) * 100L / T, σ)
    vsp = filt(get_field(TKE, "VSP", i) * 1e8L^2 / T^3, σ)
    wh = transpose(filt(get_field(a->a[:, :, z_ind], OUTPUT, "w", i) * 100L / T, σh))
    
    b = filt(get_field(a->a[:, 1, :], OUTPUT, "b", i), σ)
    b_dfm = filt(get_field(DFM, "b_dfm", i), σ)
    bh = transpose(filt(get_field(a->a[:, :, z_ind], OUTPUT, "b", i), σh))
    
    ax_w_kw = (;
        xlabel=L"x / \text{km}",
        ylabel=L"z / \text{m}", 
        limits=(-2.2sp.L * L/1000, -0sp.L * L/1000, -1.2sp.H * L, 0*L),
        xlabelsize=16,
        ylabelsize=16,
        yticks=WilkinsonTicks(3; k_max=5),
        ax_w_kw...
    )
    
    ax_vsp_kw = (;
        xlabel=L"x / \text{km}",
        ylabel=L"z / \text{m}", 
        limits=(-2.2sp.L * L/1000, -0sp.L * L/1000, -1.2sp.H * L, 0*L),
        xlabelsize=16,
        ylabelsize=16,
        yticks=WilkinsonTicks(3; k_max=5),
        ax_vsp_kw...
    )
    
    ax_wh_kw = (;
        xlabel=L"y / \text{km}", 
        ylabel=L"x / \text{km}", 
        limits=(-5sp.L * L/1000, 5sp.L * L/1000, -2.2sp.L * L/1000, -0sp.L * L/1000), 
        xlabelsize=16, 
        ylabelsize=16,
        ax_wh_kw...
    )
    
    ax_w = Axis(fig[1, 1]; ax_w_kw...)
    ax_vsp = Axis(fig[1, 3]; ax_vsp_kw...)
    ax_wh = Axis(fig[2, 1:3]; ax_wh_kw...)
    
    hideydecorations!.(ax_vsp; ticks=false)
    
    w_max = maximum(abs, wh)
    vsp_max = maximum(abs, vsp)
    
    ht_w_kw = (;
        colormap=w_colors,
        colorrange=(-w_max, w_max),
        lowclip=w_colors[1],
        highclip=w_colors[end],
        ht_w_kw...
    )
    
    ht_vsp_kw = (;
        colormap=vsp_colors,
        colorrange=(-vsp_max, vsp_max),
        lowclip=vsp_colors[1],
        highclip=vsp_colors[end],
        ht_vsp_kw...
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
    
    ht_w = heatmap!(ax_w, xsᶜ*L/1000, zsᶠ*L, w; ht_w_kw...)
    ht_vsp = heatmap!(ax_vsp, xsᶜ*L/1000, zsᶠ*L, vsp; ht_vsp_kw...)
    ht_wh = heatmap!(ax_wh, ysᶜ*L/1000, xsᶜ*L/1000, wh; ht_wh_kw...)
    
    contour!(ax_w, xsᶜ*L/1000, zsᶜ*L, b; ct_b_kw...)
    contour!(ax_vsp, xsᶜ*L/1000, zsᶜ*L, b_dfm; ct_b_kw...)
    contour!(ax_wh, ysᶜ*L/1000, xsᶜ*L/1000, bh; ct_b_kw...)
    
    lines!(ax_w, [ax_w_kw.limits[1], ax_w_kw.limits[2]], [h, h]*L; color=(:red, 0.5), linestyle=:dash)
    
    Colorbar(fig[1, 2], ht_w; label=L"w / \text{cms}^{-1}")
    Colorbar(fig[1, 4], ht_vsp; label=L"-10^8\overline{\mathbf{u}}_z \cdot \overline{w'\mathbf{u}'} / \text{m}^{2}\text{s}^{-3}")
    Colorbar(fig[2, 4], ht_wh; label=L"w / \text{cms}^{-1}")
    
    #colgap!(fig.layout, 1, 30)
    #colgap!(fig.layout, 2, 30)
    
    subfig_label!(fig[1, 1], 1)
    subfig_label!(fig[1, 3], 2)
    subfig_label!(fig[2, 1:3], 3)
    
    fig
end