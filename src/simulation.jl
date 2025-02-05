#= 
simulation.jl
    Create a simulation
    Run using
        julia ... -- simulation.jl output_folder stop_time save_interval Ro Ri α β Nx Ny Nz Q ampl init_time init_rate
=#

ENV["JULIA_SCRATCH_TRACK_ACCESS"] = 0

using Oceananigans
using JLD2


output_folder = ARGS[1]

stop_time = parse(Float64, ARGS[2])
save_interval = parse(Float64, ARGS[3])
input_parameters = let
    Ro = parse(Float64, ARGS[4])
    Ri = parse(Float64, ARGS[5])
    α = parse(Float64, ARGS[6])
    β = parse(Float64, ARGS[7])

    Nx = parse(Int, ARGS[8])
    Ny = parse(Int, ARGS[9])
    Nz = parse(Int, ARGS[10])
    
    Q = parse(Float64, ARGS[11])
    noise_amplitude = parse(Float64, ARGS[12])
    init_time = parse(Float64, ARGS[13])
    init_rate = parse(Float64, ARGS[14])
    
    preinit = parse(Bool, ARGS[15])
    preinit_path = length(ARGS) > 15 ? ARGS[16] : nothing
    
    (; Ro, Ri, α, β, Nx, Ny, Nz, Q, noise_amplitude, init_time, init_rate, preinit, preinit_path)
end
# Build parameters
include("parameters.jl")
const sp = create_simulation_parameters(input_parameters)
@info sp

# Filament state is defined here
include("filament_state.jl")

include("preinit_state.jl")

# Variable z spacing
include("z_faces.jl")

grid = RectilinearGrid(GPU(),
    size=(sp.Nx, sp.Ny, sp.Nz),
    x=(-sp.Lx/2, sp.Lx/2),
    y=(-sp.Ly/2, sp.Ly/2),
    z=z_faces(sp),
    topology=(Periodic, Periodic, Bounded),
    halo=(5, 5, 5)
)

@info grid

tracers = (:b, )
include("boundary_conditions.jl")
include("forcing.jl")

@info "Creating model"
model = NonhydrostaticModel(;
    grid,
    coriolis = FPlane(sp.f),
    clock=Clock(time=-sp.init_time),
    closure=SmagorinskyLilly(; Pr=1),
#    advection=WENO(; order=9),
    buoyancy=BuoyancyTracer(),
    tracers,
    boundary_conditions,
    forcing,
    hydrostatic_pressure_anomaly=CenterField(grid)
)

@info model

@info "Setting model state"

if sp.preinit_path == nothing
    set!(model; u=u₀, v=v₀, w=w₀, b=b₀)
else
    @info "Reading pre-init state from file"
    preinit_state = get_preinit_state(sp)
    set!(model; preinit_state...)
end

# Create a default output that saves the average state of the simulation u, v, w, b fields
@info "Creating simulation"
simulation = Simulation(model; Δt=1e-4/sp.f, stop_time)

# Progress and timestepper wizards
progress(sim) = @info string("Iteration: ", iteration(sim), ", time: ", time(sim), ", step: ", sim.Δt)
wizard = TimeStepWizard(; cfl=0.5, max_Δt=0.01/sp.f)
simulation.callbacks[:progress] = Callback(progress, TimeInterval((stop_time + sp.init_time) / 100))
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

if !isdir("$output_folder")
    mkdir(output_folder)
end
@info "Creating output writers at $output_folder"

# Fields to save
(u, v, w) = model.velocities
b = model.tracers.b
p = model.pressures.pHY′ + model.pressures.pNHS
ν = model.diffusivity_fields.νₑ

simulation.output_writers[:output] = JLD2OutputWriter(
    model,
    (; u, v, w, b, p, ν),
    filename = joinpath(output_folder, "output.jld2"),
    schedule = TimeInterval(save_interval),
    overwrite_existing = false,
    with_halos=true
)

simulation.output_writers[:checkpointer] = Checkpointer(model;
    schedule=TimeInterval(50save_interval),
    dir=output_folder,
    cleanup=true,
    overwrite_existing=false,
    verbose=true
)

# Save parameters to a file
jldopen(joinpath(output_folder, "parameters.jld2"), "w") do file
    file["simulation"] = sp
end

@info simulation
run!(simulation; pickup=true)
