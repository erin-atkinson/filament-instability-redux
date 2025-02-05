
function psi_balance_figure(
        foldernames,
        titles;
        fig_kw=(; ),
        ax_kw=(; ),
        ln_c_kw=(; ),
        ln_l_kw=(; ),
        ln_g_kw=(; ),
        ln_lg_kw=(; ),
        σ=0
    )
    
    L, T = scales()
    psi_scale = 1e6L^2 / T^3
    
    fig_kw = (; 
        size=(750, 250), 
        fig_kw...
    )
    fig = Figure(; fig_kw...)
    
    ax_kw = (;
        limits=(0, 4, -5, 5),
        xlabel=L"ft / 2\pi", 
        ylabel=L"10^6A /\text{m}^2\text{s}^{-3}",
        xlabelsize=16,
        ylabelsize=14, 
        titlesize=16,
        ax_kw...
    )
    
    ln_c_kw = (;
        color=(:black, 0.5),
        ln_c_kw...
    )
    
    map(enumerate(foldernames)) do (i, foldername)
        title = titles[i]
        ax = Axis(fig[1, i]; ax_kw..., title)
        
        # get data
        ts, c, l, g, lg = psi_data(foldername, σ)
        
        # plot the lines 
        ln_c = lines!(ax, ts, c .* psi_scale; ln_c_kw...)
        ln_l = lines!(ax, ts, l .* psi_scale; ln_l_kw...)
        ln_g = lines!(ax, ts, g .* psi_scale; ln_g_kw...)
        ln_lg = lines!(ax, ts, lg .* psi_scale; ln_lg_kw...)
        
        lns = [ln_c, ln_l, ln_g, ln_lg]
        labels = [L"{\text{d}^2C_{\text{IML}}/{\text{d}t^2}}", L"L", L"G", L"L + G"]
        
        
        subfig_label!(fig[1, i], i)
        i == length(foldernames) && Legend(fig[1, i+1], lns, labels, L"A")
        i != 1 && hideydecorations!(ax, grid=false, ticks=false)
    end
    
    fig
end

function psi_data(foldername, σ)
    sp = simulation_parameters(foldername)
    iterations, ts = iterations_times(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    
    domain = (-sp.Lx, 0, -0.85sp.H, -0.15sp.H)
    x_indices = findall(domain[1] .< xsᶜ .< domain[2])
    z_indices = findall(domain[3] .< zsᶜ .< domain[4])
    # Compute the integrals

    zsᶠ = reshape(zsᶠ, 1, length(zsᶠ))
    ml_int(field) = sum(field[x_indices, z_indices] .* (zsᶠ[:, z_indices.+1] .- zsᶠ[:, z_indices]) .* sp.Lx * sp.Ly / sp.Nx)
    
    PSI = joinpath(foldername, "PSI.jld2")
    
    Lψ = filt(timeseries_of(ml_int, PSI, "Lψ", iterations), σ)
    ∇²ψ_tt = filt(timeseries_of(ml_int, PSI, "∇²ψ_tt", iterations), σ)
    Gh = filt(timeseries_of(ml_int, PSI, "Gh", iterations), σ)
    Gv = filt(timeseries_of(ml_int, PSI, "Gv", iterations), σ)
    
    sp.f * ts / 2π, -∇²ψ_tt, Lψ, Gh .+ Gv, Lψ .+ Gh .+ Gv
end