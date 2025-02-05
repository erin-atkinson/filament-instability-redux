using JLD2
using ImageFiltering: imfilter, Kernel.gaussian
using OffsetArrays: no_offset_view
using CairoMakie
using CairoMakie.Colors
# -------------------------------------------------------------


# -------------------------------------------------------------
nov = no_offset_view

# For changing length and time scales
@inline scales() = (; L=1000, T=10000)

# Histogram
# Is there not an existing function for this?
@inline function bin_counts(data, bins)
    map(bins[1:end-1], bins[2:end]) do l, r
        count(x->l<=x<r, data)
    end
end

# Shortcut for filtering...
@inline filt(a, σ::Tuple) = imfilter(a, gaussian(σ))
@inline filt(a, σ...) = filt(a, σ)

# If only one σ, repeat for length of a
@inline filt(a::A, σ::Integer) where {T, n, A<:Array{T, n}} = filt(a, (repeat([σ], n)..., ))

@inline halos(file) = file["grid/Hx"], file["grid/Hz"], file["grid/Hz"]

# Take a vector of nd arrays and turn it into a single n+1d array
function expand(v::V) where {T, n, A<:Array{T, n}, V<:Vector{A}}
    l = length(v)
    s = size(v[1])
    new_size = (l, s...)
    a = zeros(T, new_size)
    
    inds = [Colon() for x in s] # There must be some syntax I'm missing
    for i in 1:l
        a[i, inds...] .= v[i]
    end
    a
end

expand(v::V) where {T, V<:Vector{T}} = v
# -------------------------------------------------------------


# -------------------------------------------------------------
function subfig_label!(scene, label; fontsize=14)
    gl = GridLayout(scene, 
        tellwidth = false, 
        tellheight = false, 
        halign = :left, 
        valign = :top,
    )
    Box(gl[1, 1], color = :white, strokecolor = :black, strokewidth = 0)
    Label(gl[1, 1], label;
        fontsize,
        padding = (1, 1, 2, 2),
    )
end

# Pass an integer to get that letter 
function subfig_label!(scene, label::Integer; fontsize=14)
    char = 'a' + (label - 1)
    subfig_label!(scene, "$char)"; fontsize=14)
end
# -------------------------------------------------------------


# -------------------------------------------------------------
# Grid, time and parameter retrieval
@inline function grid_nodes(file; halo=false)
    xsᶜ = nov(file["grid/xᶜᵃᵃ"])
    xsᶠ = nov(file["grid/xᶠᵃᵃ"])

    ysᶜ = nov(file["grid/yᵃᶜᵃ"])
    ysᶠ = nov(file["grid/yᵃᶠᵃ"])

    zsᶜ = nov(file["grid/zᵃᵃᶜ"])
    zsᶠ = nov(file["grid/zᵃᵃᶠ"])

    if !halo
        Hx, Hy, Hz = halos(file)

        xsᶜ = xsᶜ[(Hx+1):(end-Hx)]
        xsᶠ = xsᶠ[(Hx+1):(end-Hx)]
        ysᶜ = ysᶜ[(Hy+1):(end-Hy)]
        ysᶠ = ysᶠ[(Hy+1):(end-Hy)]
        zsᶜ = zsᶜ[(Hz+1):(end-Hz)]
        zsᶠ = zsᶠ[(Hz+1):(end-Hz)]
    end

    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ
end

@inline function grid_nodes(foldername::String; halo=false)
    
    filename = joinpath(foldername, "output.jld2")
    !isfile(filename) && (filename = joinpath(foldername, "DFM.jld2"))
    
    jldopen(filename) do file
        grid_nodes(file; halo)
    end 
end

@inline function iterations_times(file)
    
    iterations = keys(file["timeseries/t"])
    ts = map(iteration->file["timeseries/t/$iteration"], iterations)
    
    iterations, ts
end

@inline function iterations_times(foldername::String)
    
    filename = joinpath(foldername, "output.jld2")
    !isfile(filename) && (filename = joinpath(foldername, "DFM.jld2"))
    
    jldopen(iterations_times, filename)
end

@inline function simulation_parameters(foldername)
    
    filename = joinpath(foldername, "parameters.jld2")
    
    jldopen(filename) do file
        file["simulation"]
    end
end
# -------------------------------------------------------------


# -------------------------------------------------------------
# Get field from the JLD2 file
@inline function get_field(file, field, iteration; halo=false)
    
    jldpath = "timeseries/$field/$iteration"
    data = nov(file[jldpath])
    
    # Remove singleton dimensions and halo points
    Hs = halo ? zeros(3) : halos(file)
    indices = map(Hs, size(data)) do H, s
        s==1 ? 1 : (1+H):(s-H)
    end
    
    data[indices...]
end

# Get from a string
@inline function get_field(filename::String, field, iteration)
    jldopen(file->get_field(file, field, iteration), filename)
end

# Apply a function to the field after getting
@inline function get_field(f::Function, args...)
    f(get_field(args...))
end
# -------------------------------------------------------------


# -------------------------------------------------------------
# Reduce Nd field using func (Nd->Md), then find it for all times
@inline function timeseries_of(func, file, field, iterations)
    data = map(iteration->get_field(func, file, field, iteration), iterations)
    expand(data)
end

# From a string
@inline function timeseries_of(func, filename::String, args...)
    jldopen(filename) do file
        timeseries_of(func, file, args...)
    end
end
# -------------------------------------------------------------


# -------------------------------------------------------------
# Time average of a function
@inline function time_average_of(func, file, field, iterations)
    mapreduce(.+, iterations) do iteration
        get_field(func, file, field, iteration)
    end ./ length(iterations)
end

# From a string
@inline function time_average_of(func, filename::String, args...)
    jldopen(filename) do file
        time_average_of(func, file, args...)
    end
end
# -------------------------------------------------------------


# -------------------------------------------------------------
for i in [:x, :y, :z], j in [:ᶜ, :ᶠ]
    fname = Symbol(i, j, :bounds)
    nodes = Symbol(i, :s, j)
    @eval function $fname(foldername, args::Tuple)
        xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
        map(arg->argmin(abs.($nodes .- arg)), args)
    end
    @eval $fname(foldername, args...) = $fname(foldername, args)
    @eval $fname(foldername, i) = $fname(foldername, (i, ))[1]
end

function tbounds(foldername, args::Tuple)
    iterations, times = iterations_times(foldername)
    map(arg->argmin(abs.(times .- arg)), args)
end
tbounds(foldername, args...) = tbounds(foldername, args)
tbounds(foldername, i) = tbounds(foldername, (i, ))[1]
# -------------------------------------------------------------
