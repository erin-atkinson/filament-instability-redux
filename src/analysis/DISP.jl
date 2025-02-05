using Oceananigans.Operators

# Frontogenetic tendency equation
include("terms/terms.jl")
include("terms/TKE_terms.jl")



b_dfm = dfm(b)

@inline function b_x_func(i, j, k, grid, b_dfm)
    return ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b_dfm)
end

DISP_op = KernelFunctionOperation{Center, Center, Center}(DISP_func, grid, u, v, w, ν)
b_x_op = KernelFunctionOperation{Center, Nothing, Center}(b_x_func, grid, b_dfm)

DISP = Field(DISP_op)
DISP_dfm = dfm(DISP)

b_x = Field(b_x_op)

outputs = (; b_x, DISP_dfm)

function update_outputs!(outputs)
    compute!(b_dfm)
    compute!(DISP)
    map(compute!, outputs)
end

macro postpostprocess()
    
end