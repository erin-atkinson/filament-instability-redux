function wVSP_video(
        foldername,
        output_filename,
        h,
        frames=Colon(); 
        fig_kw=(; ), 
        ax_w_kw=(; title=""),
        ax_vsp_kw=(; title=""),
        ax_wh_kw=(; title=""),
        ht_w_kw=(; ),
        ht_vsp_kw=(; ),
        ht_wh_kw=(; ),
        ct_b_kw=(; ),
        σ=0,
        σh=0
    )
    
    L, T = scales()
    b_levels = range(-300, -227, 16)
    w_colormap = to_colormap(:balance)
    vsp_colormap = reverse(to_colormap(:curl))
    
    fig_kw = (;
        size=(750, 500),
        fig_kw...
    )
    
    fig = Figure(; fig_kw...)
    
    sp = simulation_parameters(foldername)
    iterations, ts = iterations_times(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    
    z_ind = zᶠbounds(foldername, h)
    
    DFM_filename = joinpath(foldername, "DFM.jld2")
    TKE_filename = joinpath(foldername, "TKE.jld2")
    OUTPUT_filename = joinpath(foldername, "output.jld2")
    # Find maximum for a constant colorrange
    vsp_max = maximum(timeseries_of(a->maximum(abs, filt(a, σ)) * 1e8L^2/T^3, TKE_filename, "VSP", iterations[frames]))
    w_max = maximum(timeseries_of(a->maximum(abs, filt(a[:, 1, :], σ)) * 1e2L/T, OUTPUT_filename, "w", iterations[frames]))
    wh_max = maximum(timeseries_of(a->maximum(abs, filt(a[:, :, z_ind], σh)) * 1e2L/T, OUTPUT_filename, "w", iterations[frames]))
    
    vsp_clip = :colorrange in keys(ht_vsp_kw) && ht_vsp_kw.colorrange[1] < vsp_max
    w_clip = :colorrange in keys(ht_w_kw) && ht_w_kw.colorrange[1] < w_max
    wh_clip = :colorrange in keys(ht_wh_kw) && ht_wh_kw.colorrange[1] < wh_max
    
    frames = let a = 1:length(iterations)
        a[frames]
    end
    frame = Observable(frames[1])
    iteration = @lift iterations[$frame]
    t = @lift ts[$frame]
    
    # Set each title
    d_vsp_title = :title in keys(ax_vsp_kw) ? ax_vsp_kw.title : ""
    vsp_title = @lift let t_val = rpad(round(sp.f * $t / 2π; digits=2), 4, '0')
        hr_val = rpad(round(T * $t / 3600; digits=2), 5, '0')
        L"%$d_vsp_title $ft / 2\pi = %$t_val \quad t = %$hr_val~\text{hr}$"
    end
    d_w_title = :title in keys(ax_w_kw) ? ax_w_kw.title : ""
    w_title = @lift let t_val = rpad(round(sp.f * $t / 2π; digits=2), 4, '0')
        hr_val = rpad(round(T * $t / 3600; digits=2), 5, '0')
        L"%$d_w_title $ft / 2\pi = %$t_val \quad t = %$hr_val~\text{hr}$"
    end
    d_wh_title = :title in keys(ax_wh_kw) ? ax_wh_kw.title : ""
    wh_title = @lift let t_val = rpad(round(sp.f * $t / 2π; digits=2), 4, '0')
        hr_val = rpad(round(T * $t / 3600; digits=2), 5, '0')
        L"%$d_wh_title $ft / 2\pi = %$t_val \quad t = %$hr_val~\text{hr}$"
    end
    
    ax_vsp_kw = (;
        xlabel=L"x/\text{km}",
        ylabel=L"z/\text{m}",
        limits=(-2.2sp.L * L/1000, -0sp.L * L/1000, -1.2sp.H * L, 0*L),
        ax_vsp_kw...,
        title=""
    )
    ax_w_kw = (; 
        xlabel=L"x/\text{km}",
        ylabel=L"z/\text{m}",
        limits=(-2.2sp.L * L/1000, -0sp.L * L/1000, -1.2sp.H * L, 0*L),
        ax_w_kw...,
        title=w_title
    )
    ax_wh_kw = (; 
        xlabel=L"y/\text{km}",
        ylabel=L"x/\text{km}",
        limits=(-5sp.L * L/1000, 5sp.L * L/1000, -2.2sp.L * L/1000, -0sp.L * L/1000),
        ax_wh_kw...,
        title=""
    )
    
    ht_vsp_kw=(;
        colorrange=(-vsp_max, vsp_max),
        colormap=vsp_colormap,
        highclip= vsp_clip ? vsp_colormap[end] : nothing,
        lowclip= vsp_clip ? vsp_colormap[1] : nothing,
        ht_vsp_kw...
    )
    ht_w_kw=(;
        colorrange=(-w_max, w_max),
        colormap=w_colormap,
        highclip= w_clip ? w_colormap[end] : nothing,
        lowclip= w_clip ? w_colormap[1] : nothing,
        ht_w_kw...
    )
    ht_wh_kw=(;
        colorrange=(-wh_max, wh_max),
        colormap=w_colormap,
        highclip= wh_clip ? w_colormap[end] : nothing,
        lowclip= wh_clip ? w_colormap[1] : nothing,
        ht_wh_kw...
    )
    ct_b_kw=(;
        levels=b_levels,
        color=(:black, 0.5),
        linewidth=2,
        ct_b_kw...
    )
    
    ax_vsp = Axis(fig[1, 3]; ax_vsp_kw...)
    ax_w = Axis(fig[1, 1]; ax_w_kw...)
    ax_wh = Axis(fig[2, 1:3]; ax_wh_kw...)
    hideydecorations!.(ax_vsp; ticks=false)
    
    DFM = jldopen(DFM_filename)
    TKE = jldopen(TKE_filename)
    OUTPUT = jldopen(OUTPUT_filename)
    
    vsp = @lift get_field(a->filt(a, σ) * 1e8L^2/T^3, TKE, "VSP", $iteration)
    w = @lift get_field(a->filt(a[:, 1, :], σ) * 1e2L/T, OUTPUT, "w", $iteration)
    wh = @lift get_field(a->transpose(filt(a[:, :, z_ind], σh)) * 1e2L/T, OUTPUT, "w", $iteration)
    
    b_dfm = @lift get_field(a->filt(a, σ), DFM, "b_dfm", $iteration)
    b = @lift get_field(a->filt(a[:, 1, :], σ), OUTPUT, "b", $iteration)
    bh = @lift get_field(a->transpose(filt(a[:, :, z_ind], σh)), OUTPUT, "b", $iteration)
    
    ht_vsp = heatmap!(ax_vsp, xsᶜ*L/1000, zsᶜ*L, vsp; ht_vsp_kw...)
    ht_w = heatmap!(ax_w, xsᶜ*L/1000, zsᶠ*L, w; ht_w_kw...)
    ht_wh = heatmap!(ax_wh, ysᶜ*L/1000, xsᶜ*L/1000, wh; ht_wh_kw...)
    
    contour!(ax_vsp, xsᶜ*L/1000, zsᶜ*L, b_dfm; ct_b_kw...)
    contour!(ax_w, xsᶜ*L/1000, zsᶜ*L, b; ct_b_kw...)
    contour!(ax_wh, ysᶜ*L/1000, xsᶜ*L/1000, bh; ct_b_kw...)
    
    Colorbar(fig[1, 4], ht_vsp; label=L"-10^8\overline{\mathbf{u}}_z \cdot \overline{w'\mathbf{u}'} / \text{m}^{2}\text{s}^{-3}")
    Colorbar(fig[1, 2], ht_w; label=L"w(x, 0, z) / \text{cms}^{-1}")
    Colorbar(fig[2, 4], ht_wh; label=L"w(x, y, z=%$(round(h * L; sigdigits=2))~\text{m}) / \text{cms}^{-1}")
    
    subfig_label!(fig[1, 1], 1)
    subfig_label!(fig[1, 3], 2)
    subfig_label!(fig[2, 1:3], 3)
    
    record(fig, output_filename, frames; framerate=12) do i
        frame[] = i
        print("$output_filename: $(frames[1])->$i->$(frames[end])\r")
    end
    println("")
    
    close(DFM)
    close(TKE)
    close(OUTPUT)
    
    fig
end