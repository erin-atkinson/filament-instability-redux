# postinit_state.jl

using OffsetArrays: no_offset_view
using Statistics: mean

function get_preinit_state(sp; n=0)
    
    filename = joinpath(sp.preinit_path, "output.jld2")
    
    iterations = jldopen(file->keys(file["timeseries/t"]), filename) 
    grid = jldopen(file->file["serialized/grid"], filename) 
    
    n = n == 0 ? length(iterations) : n
    iteration = iterations[n]
    
    u = interior(FieldTimeSeries(filename, "u"; iterations=[iteration])[1])
    v = interior(FieldTimeSeries(filename, "v"; iterations=[iteration])[1])
    w = interior(FieldTimeSeries(filename, "w"; iterations=[iteration])[1])
    b = interior(FieldTimeSeries(filename, "b"; iterations=[iteration])[1])
    
    u = u .- mean(u; dims=2)
    v = v .- mean(v; dims=2)
    w = w .- mean(w; dims=2)
    b = b .- mean(b; dims=2)

    
    xsᶜ = no_offset_view(grid.xᶜᵃᵃ)[6:end-5]
    xsᶠ = no_offset_view(grid.xᶠᵃᵃ)[6:end-5]
    ysᶜ = no_offset_view(grid.yᵃᶜᵃ)[6:end-5]
    ysᶠ = no_offset_view(grid.yᵃᶠᵃ)[6:end-5]
    zsᶜ = no_offset_view(grid.zᵃᵃᶜ)[6:end-5]
    zsᶠ = no_offset_view(grid.zᵃᵃᶠ)[6:end-5]
    
    # Scale?
    u_init = sp.noise_amplitude .* u
    w_init = sp.noise_amplitude .* w
    
    # Add on the filament
    v_filament = [v₀(x, y, z) for x in xsᶜ, y in ysᶠ, z in zsᶜ] .+ (sp.noise_amplitude .* v)
    b_filament = [b₀(x, y, z) for x in xsᶜ, y in ysᶜ, z in zsᶜ] .+ (sp.noise_amplitude .* b)
    
    return (; u=u_init, w=w_init, v=v_filament, b=b_filament)
end
