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
        size=(500, 1000),
        fig_kw...
    )
    
    fig = Figure(; fig_kw...)
    
    sp = simulation_parameters(foldername)
    iterations, ts = iterations_times(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    
    filename = joinpath(foldername, "DFM.jl")
    # Find maximum for a constant colorrange
    u_max = timeseries_of(a->maximum(abs, filt(a, σ)) * L/T, filename, "u_dfm", iterations[frames])
    v_max = timeseries_of(a->maximum(abs, filt(a, σ)) * L/T, filename, "v_dfm", iterations[frames])
    w_max = timeseries_of(a->maximum(abs, filt(a, σ)) * L/T, filename, "w_dfm", iterations[frames])
    
    u_clip = ht_w_kw.colorrange[1] < u_max
    v_clip = ht_v_kw.colorrange[1] < v_max
    w_clip = ht_w_kw.colorrange[1] < w_max
    
    frames = let a = 1:length(iterations)
        a[frames]
    end
    frame = Observable(frames[1])
    iteration = @lift iterations[$frame]
    t = @lift ts[$frame]
    
    # Set each title
    u_title = @lift let t_val = round(sp.f * $t / 2π; digits=2)
        L"%$(ax_u_kw.title) $ft / 2\pi = %$t_val$"
    end
    v_title = @lift let t_val = round(sp.f * $t / 2π; digits=2)
        L"%$(ax_v_kw.title) $ft / 2\pi = %$t_val$"
    end
    w_title = @lift let t_val = round(sp.f * $t / 2π; digits=2)
        L"%$(ax_w_kw.title) $ft / 2\pi = %$t_val$"
    end
    
    ax_u_kw = (;
        xlabel=L"x/\text{km}",
        ylabel=L"z/\text{m}",
        ax_u_kw...,
        title=u_title
    )
    ax_v_kw = (; 
        xlabel=L"x/\text{km}",
        ylabel=L"z/\text{m}",
        ax_v_kw...,
        title=v_title
    )
    ax_w_kw = (; 
        xlabel=L"x/\text{km}",
        ylabel=L"z/\text{m}",
        ax_w_kw...,
        title=w_title
    )
    
    ht_u_kw=(;
        colorrange=(-u_max, u_max),
        colormap,
        ht_u_kw...
    )
    ht_v_kw=(;
        colorrange=(-v_max, v_max),
        colormap,
        ht_v_kw...
    )
    ht_w_kw=(;
        colorrange=(-w_max, w_max),
        colormap,
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
    
    u = @lift get_field(a->filt(a, σ), DFM, "u_dfm", $iteration)
    v = @lift get_field(a->filt(a, σ), DFM, "v_dfm", $iteration)
    w = @lift get_field(a->filt(a, σ), DFM, "w_dfm", $iteration)
    b = @lift get_field(a->filt(a, σ), DFM, "b_dfm", $iteration)
    
    ht_u = heatmap!(ax_u, xsᶠ, zsᶜ, u; ht_u_kw)
    ht_v = heatmap!(ax_v, xsᶜ, zsᶜ, v; ht_v_kw)
    ht_w = heatmap!(ax_w, xsᶜ, zsᶠ, w; ht_w_kw)
    
    contour!(ax_u, xsᶜ, zsᶜ, b; ct_b_kw)
    contour!(ax_v, xsᶜ, zsᶜ, b; ct_b_kw)
    contour!(ax_w, xsᶜ, zsᶜ, b; ct_b_kw)
    
    Colorbar(fig[1, 2], ht_u; label=L"u / \text{ms}^{-1}")
    Colorbar(fig[2, 2], ht_v; label=L"v / \text{ms}^{-1}")
    Colorbar(fig[3, 2], ht_w; label=L"w / \text{ms}^{-1}")
    
    record(fig, output_filename, frames) do i
        frame[] = i
        print("$(frames[1])->$i->$(frames[end])")
    end
    
    close(DFM)
    
    fig
end