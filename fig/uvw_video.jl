function uvw_video(
        foldername,
        output_filename,
        frames=Colon(); 
        fig_kw=(; ), 
        ax_u_kw=(; title=""),
        ax_v_kw=(; title=""),
        ax_w_kw=(; title=""),
        ht_u_kw=(; ),
        ht_v_kw=(; ),
        ht_w_kw=(; ),
        ct_b_kw=(; ),
        σ=0
    )
    
    L, T = scales()
    b_levels = range(-300, -227, 16)
    colormap = to_colormap(:balance)
    
    fig_kw = (;
        size=(500, 750),
        fig_kw...
    )
    
    fig = Figure(; fig_kw...)
    
    sp = simulation_parameters(foldername)
    iterations, ts = iterations_times(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    
    filename = joinpath(foldername, "DFM.jld2")
    # Find maximum for a constant colorrange
    u_max = maximum(timeseries_of(a->maximum(abs, filt(a, σ)) * L/T, filename, "u_dfm", iterations[frames]))
    v_max = maximum(timeseries_of(a->maximum(abs, filt(a, σ)) * L/T, filename, "v_dfm", iterations[frames]))
    w_max = maximum(timeseries_of(a->maximum(abs, filt(a, σ)) * L/T, filename, "w_dfm", iterations[frames]))
    
    u_clip = :colorrange in keys(ht_u_kw) && ht_u_kw.colorrange[1] < u_max
    v_clip = :colorrange in keys(ht_v_kw) && ht_v_kw.colorrange[1] < v_max
    w_clip = :colorrange in keys(ht_w_kw) && ht_w_kw.colorrange[1] < w_max
    
    frames = let a = 1:length(iterations)
        a[frames]
    end
    frame = Observable(frames[1])
    iteration = @lift iterations[$frame]
    t = @lift ts[$frame]
    
    # Set each title
    d_u_title = :title in keys(ax_u_kw) ? ax_u_kw.title : ""
    u_title = @lift let t_val = round(sp.f * $t / 2π; digits=2)
        L"%$d_u_title $ft / 2\pi = %$t_val$"
    end
    d_v_title = :title in keys(ax_v_kw) ? ax_v_kw.title : ""
    v_title = @lift let t_val = round(sp.f * $t / 2π; digits=2)
        L"%$d_v_title $ft / 2\pi = %$t_val$"
    end
    d_w_title = :title in keys(ax_w_kw) ? ax_w_kw.title : ""
    w_title = @lift let t_val = round(sp.f * $t / 2π; digits=2)
        L"%$d_w_title $ft / 2\pi = %$t_val$"
    end
    
    ax_u_kw = (;
        xlabel=L"x/\text{km}",
        ylabel=L"z/\text{m}",
        limits=(-3*L*sp.L/1000, 3*L*sp.L/1000, -1.2sp.H * L, 0),
        ax_u_kw...,
        title=u_title
    )
    ax_v_kw = (; 
        xlabel=L"x/\text{km}",
        ylabel=L"z/\text{m}",
        limits=(-3*L*sp.L/1000, 3*L*sp.L/1000, -1.2sp.H * L, 0),
        ax_v_kw...,
        title=v_title
    )
    ax_w_kw = (; 
        xlabel=L"x/\text{km}",
        ylabel=L"z/\text{m}",
        limits=(-3*L*sp.L/1000, 3*L*sp.L/1000, -1.2sp.H * L, 0),
        ax_w_kw...,
        title=w_title
    )
    
    ht_u_kw=(;
        colorrange=(-u_max, u_max),
        colormap,
        highclip= u_clip ? colormap[end] : nothing,
        lowclip= u_clip ? colormap[1] : nothing,
        ht_u_kw...
    )
    ht_v_kw=(;
        colorrange=(-v_max, v_max),
        colormap,
        highclip= v_clip ? colormap[end] : nothing,
        lowclip= v_clip ? colormap[1] : nothing,
        ht_v_kw...
    )
    ht_w_kw=(;
        colorrange=(-w_max, w_max),
        colormap,
        highclip= w_clip ? colormap[end] : nothing,
        lowclip= w_clip ? colormap[1] : nothing,
        ht_w_kw...
    )
    ct_b_kw=(;
        levels=b_levels,
        color=(:black, 0.5),
        linewidth=2,
        ct_b_kw...
    )
    
    ax_u = Axis(fig[1, 1]; ax_u_kw...)
    ax_v = Axis(fig[2, 1]; ax_v_kw...)
    ax_w = Axis(fig[3, 1]; ax_w_kw...)
    
    DFM = jldopen(filename)
    
    u = @lift get_field(a->filt(a, σ) * L/T, DFM, "u_dfm", $iteration)
    v = @lift get_field(a->filt(a, σ) * L/T, DFM, "v_dfm", $iteration)
    w = @lift get_field(a->filt(a, σ) * L/T, DFM, "w_dfm", $iteration)
    b = @lift get_field(a->filt(a, σ), DFM, "b_dfm", $iteration)
    
    ht_u = heatmap!(ax_u, xsᶠ*L/1000, zsᶜ*L, u; ht_u_kw...)
    ht_v = heatmap!(ax_v, xsᶜ*L/1000, zsᶜ*L, v; ht_v_kw...)
    ht_w = heatmap!(ax_w, xsᶜ*L/1000, zsᶠ*L, w; ht_w_kw...)
    
    contour!(ax_u, xsᶜ*L/1000, zsᶜ*L, b; ct_b_kw...)
    contour!(ax_v, xsᶜ*L/1000, zsᶜ*L, b; ct_b_kw...)
    contour!(ax_w, xsᶜ*L/1000, zsᶜ*L, b; ct_b_kw...)
    
    Colorbar(fig[1, 2], ht_u; label=L"u / \text{ms}^{-1}")
    Colorbar(fig[2, 2], ht_v; label=L"v / \text{ms}^{-1}")
    Colorbar(fig[3, 2], ht_w; label=L"w / \text{ms}^{-1}")
    
    record(fig, output_filename, frames; framerate=12) do i
        frame[] = i
        print("$(frames[1])->$i->$(frames[end])\r")
    end
    
    close(DFM)
    
    fig
end