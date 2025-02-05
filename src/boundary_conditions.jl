#=
boundary_conditions.jl
    Boundary conditions for the buoyancy, cooling on the top and bottom stratification.
=#

using Oceananigans

b = FieldBoundaryConditions(;
    top=FluxBoundaryCondition(sp.Q),
    bottom=GradientBoundaryCondition(sp.N₀²)
)

boundary_conditions = (; b)