using Oceananigans.Operators

# Frontogenetic tendency equation
include("terms/terms.jl")
include("terms/fluctuation_terms.jl")
include("terms/front_terms.jl")

u_dfm = dfm(u)
w_dfm = dfm(w)
b_dfm = dfm(b)

fields = (; u, w, b)
mean_fields = (; u_dfm, w_dfm, b_dfm)

T_adv_op = KernelFunctionOperation{Center, Nothing, Center}(T_adv_func, grid, mean_fields)
Fh_op = KernelFunctionOperation{Center, Nothing, Center}(Fh_func, grid, fields, mean_fields)
Fv_op = KernelFunctionOperation{Center, Nothing, Center}(Fv_func, grid, fields, mean_fields)
∇b²_op = KernelFunctionOperation{Center, Nothing, Center}(∇b²_func, grid, mean_fields)

∇b² = Field(∇b²_op)

T_adv = Field(T_adv_op)
Fh = Field(Fh_op)
Fv = Field(Fv_op)

∇b²_adv_op = KernelFunctionOperation{Center, Nothing, Center}(∇b²_adv_func, grid, mean_fields, ∇b²)
∇b²_adv = Field(∇b²_adv_op)

outputs = (; T_adv, Fh, Fv, ∇b²)
temp_outputs = (; ∇b²_adv, ∇b²)

function update_outputs!(outputs)
    map(compute!, mean_fields)
    map(compute!, outputs)
    map(compute!, temp_outputs)
end

macro postpostprocess()
    Δt = ts[2] - ts[1]
    jldopen(temp_filename, "r") do temp_file
        
        @inline ∇b²(iter) = temp_file["timeseries/∇b²/$iter"]
        @inline ∇b²_adv(iter) = temp_file["timeseries/∇b²_adv/$iter"]
        
        for i in 1:length(iters)
            
            
            p_iter = iters[max(1, i-1)]
            iter = iters[i]
            n_iter = iters[min(length(iters), i+1)]
            
            
            D∇b²Dt = ∇b²_adv(iter) .+ (∇b²(n_iter) .- ∇b²(p_iter)) ./ (2Δt)
            jldopen(filename, "a") do file
                file["timeseries/D∇b²Dt/$iter"] = D∇b²Dt
            end
        end
    end
end