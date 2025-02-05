using Oceananigans.Operators

# This calculates linear terms in the streamfunction evolution equation
# (no non-linear, turbulence, or sgs)
# doesn't compute time derivatives

# Need to grab the base state
include("../filament_state.jl")

include("terms/terms.jl")
include("terms/fluctuation_terms.jl")
include("terms/psi_terms.jl")

v_field, b_field = let
    xs = xnodes(grid, Center(); with_halos=true)
    zs = znodes(grid, Center(); with_halos=true)
    v_data = [v₀(x, 0, z) for x in xs, y in 1:1, z in zs]
    b_data = [b₀(x, 0, z) for x in xs, y in 1:1, z in zs]
    Field{Center, Nothing, Center}(grid; data=v_data), Field{Center, Nothing, Center}(grid; data=b_data)
end 

u_dfm = dfm(u)
v_dfm = dfm(v)
w_dfm = dfm(w)
b_dfm = dfm(b)

ψ = Field(@at (Center, Nothing, Center) CumulativeIntegral(-u_dfm; dims=3))

fields = (; u, v, w, b)
mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm, ψ)


Lψ_op = KernelFunctionOperation{Center, Nothing, Center}(Lψ_func, grid, v_field, b_field, ψ, sp)
∇²ψ_op = KernelFunctionOperation{Center, Nothing, Center}(∇²ψ_func, grid, ψ)


Gh_t_op = KernelFunctionOperation{Center, Nothing, Center}(Gh_t_func, grid, fields, mean_fields, sp)
Gh_nt_op = KernelFunctionOperation{Center, Nothing, Center}(Gh_nt_func, grid, fields, mean_fields, sp)
Gv_t_op = KernelFunctionOperation{Center, Nothing, Center}(Gv_t_func, grid, fields, mean_fields, sp)
Gv_nt_op = KernelFunctionOperation{Center, Nothing, Center}(Gv_nt_func, grid, fields, mean_fields, sp)

Lψ = Field(Lψ_op)
∇²ψ = Field(∇²ψ_op)

Gh_t = Field(Gh_t_op)
Gh_nt = Field(Gh_nt_op)
Gv_t = Field(Gv_t_op)
Gv_nt = Field(Gv_nt_op)

outputs = (; ψ, Lψ, ∇²ψ)
temp_outputs = (; ∇²ψ, Gh_t, Gh_nt, Gv_t, Gv_nt)

function update_outputs!(outputs)
    map(compute!, mean_fields)
    map(compute!, temp_outputs)
    map(compute!, outputs)
end

# Custom script for processing-after-processing?
macro postpostprocess()
    @info "Post-post-processing for time derivative of ψ"
    # time derivative tt of ψ
    # continue at boundaries
    Δt = ts[2] - ts[1]
    jldopen(temp_filename, "r") do temp_file
        
        @inline ∇²ψ(iter) = temp_file["timeseries/∇²ψ/$iter"]
        
        @inline Gh_t(iter) = temp_file["timeseries/Gh_t/$iter"]
        @inline Gh_nt(iter) = temp_file["timeseries/Gh_nt/$iter"]
        @inline Gv_t(iter) = temp_file["timeseries/Gv_t/$iter"]
        @inline Gv_nt(iter) = temp_file["timeseries/Gv_nt/$iter"]
        
        for i in 1:length(iters)
            
            
            p_iter = iters[max(1, i-1)]
            iter = iters[i]
            n_iter = iters[min(length(iters), i+1)]
            
            ∇²ψ_tt = -(2 / Δt^2) * (∇²ψ(iter) - (∇²ψ(n_iter) + ∇²ψ(p_iter)) / 2)
            Gh = Gh_nt(iter) .+ (Gh_t(n_iter) .- Gh_t(p_iter)) ./ (2Δt)
            Gv = Gv_nt(iter) .+ (Gv_t(n_iter) .- Gv_t(p_iter)) ./ (2Δt)
            jldopen(filename, "a") do file
                file["timeseries/∇²ψ_tt/$iter"] = ∇²ψ_tt
                file["timeseries/Gh/$iter"] = Gh
                file["timeseries/Gv/$iter"] = Gv
            end
        end
    end
end