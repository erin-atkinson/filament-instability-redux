using Oceanostics
using Oceananigans.Operators

# Kinetic energy balance in the down-front mean fields
# Need to grab the base state
include("../filament_state.jl")

include("terms/terms.jl")
include("terms/fluctuation_terms.jl")
include("terms/TKE_terms.jl")

v_field, b_field = let
    xs = xnodes(grid, Center(); with_halos=true)
    zs = znodes(grid, Center(); with_halos=true)
    v_data = [v₀(x, 0, z) for x in xs, y in 1:1, z in zs]
    b_data = [b₀(x, 0, z) for x in xs, y in 1:1, z in zs]
    Field{Center, Nothing, Center}(grid; data=v_data), Field{Center, Nothing, Center}(grid; data=b_data)
end 

mean_fields = (; u_dfm=Field{Face, Nothing, Center}(grid), v_dfm=v_field, w_dfm=Field{Center, Nothing, Face}(grid), b_dfm=b_field)

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
    map(compute!, outputs)
end
