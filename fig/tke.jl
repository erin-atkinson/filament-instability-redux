# tke.jl
# A comparison of tke production terms across simulations

function tke_figure(
        foldernames,
        titles;
        fig_kw=(; ),
        ax_kw=(; ),
        ln_lsp_kw=(; ),
        ln_vsp_kw=(; ),
        ln_bflux_kw=(; ),
        tke₀=false,
        ref_foldername="",
    )
    
    # get normalisation
    ref_foldername = ref_foldername == "" ? foldernames[1] : ref_foldername
    tke_scale = maximum(abs, tke_data(ref_foldername)[3])
    
    fig_kw = (; 
        size=(750, 230), 
        fig_kw...
    )
    fig = Figure(; fig_kw...)
    
    ax_kw = (;
        limits=(0, 4, -0.2, 1.2),
        xlabel=L"ft / 2\pi", 
        ylabel="",
        xlabelsize=16,
        ylabelsize=14, 
        titlesize=16,
        ax_kw...
    )
    
    ln_lsp_kw = (;
        color=:red,
        ln_lsp_kw...
    )
    ln_vsp_kw = (;
        color=:green,
        ln_vsp_kw...
    )
    ln_bflux_kw = (;
        color=:blue,
        ln_bflux_kw...
    )
    
    map(enumerate(foldernames)) do (i, foldername)
        title = titles[i]
        ax = Axis(fig[1, i]; ax_kw..., title)
        
        # get data
        ts, lsp, vsp, bflux = tke_data(foldername, tke₀)
        
        # plot the lines 
        ln_lsp = lines!(ax, ts, lsp ./ tke_scale; ln_lsp_kw...)
        ln_bflux = lines!(ax, ts, bflux ./ tke_scale; ln_bflux_kw...)
        ln_vsp = lines!(ax, ts, vsp ./ tke_scale; ln_vsp_kw...)
        
        # Also plot a reference if zero
        if tke₀
            _, _, vsp_reg, _ = tke_data(foldername)
            ln_vsp_reg = lines!(ax, ts, vsp_reg ./ tke_scale; ln_vsp_kw..., linestyle=:dash)
            
            lns = [ln_lsp, ln_bflux, ln_vsp, ln_vsp_reg]
            labels = [L"\text{LSP}_0", L"\text{BFLUX}_0", L"\text{VSP}_0", L"\text{VSP}"]
        else
            lns = [ln_lsp, ln_bflux, ln_vsp]
            labels = [L"\text{LSP}", L"\text{BFLUX}", L"\text{VSP}"]
        end
        
        
        subfig_label!(fig[1, i], i)
        i == length(foldernames) && Legend(fig[1, i+1], lns, labels)
        i != 1 && hideydecorations!(ax, grid=false, ticks=false)
    end
    
    fig
end

function tke_data(foldername, tke₀=false)
    sp = simulation_parameters(foldername)
    iterations, ts = iterations_times(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    
    z_indices = findall(-2sp.H .<= zsᶜ .<= 0)
    zsᶠ = reshape(zsᶠ, 1, length(zsᶠ))
    ml_int(field) = sum(field[:, z_indices] .* (zsᶠ[:, z_indices.+1] .- zsᶠ[:, z_indices]) .* sp.Lx * sp.Ly / sp.Nx)
    
    TKE = tke₀ ? joinpath(foldername, "TKE-NL.jld2") : joinpath(foldername, "TKE.jld2")
    
    lsp = timeseries_of(ml_int, TKE, "LSP", iterations)
    vsp = timeseries_of(ml_int, TKE, "VSP", iterations)
    bflux = timeseries_of(ml_int, TKE, "BFLUX", iterations)
    
    sp.f * ts / 2π, lsp, vsp, bflux
end