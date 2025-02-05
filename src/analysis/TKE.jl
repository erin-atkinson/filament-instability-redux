using Oceanostics
using Oceananigans.Operators

# Kinetic energy balance in the down-front mean fields

include("terms/terms.jl")
include("terms/fluctuation_terms.jl")
include("terms/TKE_terms.jl")

u_dfm = dfm(u)
v_dfm = dfm(v)
w_dfm = dfm(w)
b_dfm = dfm(b)

mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm)

# Shear production
LSP_op = KernelFunctionOperation{Center, Nothing, Center}(LSP_func, grid, u, v, w, mean_fields)
LSP = Field(LSP_op)

VSP_op = KernelFunctionOperation{Center, Nothing, Center}(VSP_func, grid, u, v, w, mean_fields)
VSP = Field(VSP_op)

# Potential energy
BFLUX_op = KernelFunctionOperation{Center, Nothing, Center}(BFLUX_func, grid, w, b, mean_fields)
BFLUX = Field(BFLUX_op)

outputs = (; LSP, VSP, BFLUX)

function update_outputs!(outputs)
    map(compute!, mean_fields)
    map(compute!, outputs)
end
