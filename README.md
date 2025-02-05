# filament-instability-redux
Setup for simulation of a mixed layer filament using Oceananigans. Setup to easily change reference conditions, initialisation method, cooling etc. Hopefully this is useful later. Simulations and post processing done on SciNet Mist, Niagara systems respectively.

## Setup
Environment used on Mist (PowerPC) in `$SCRATCH/.julia-mist`
```
  [a9b6321e] Atomix v0.1.0 ⚲
  [052768ef] CUDA v5.3.5 ⚲
  [033835bb] JLD2 v0.5.10
  [9e8cae18] Oceananigans v0.93.3 ⚲
  [6fe1bfb0] OffsetArrays v1.15.0
  [276daf66] SpecialFunctions v2.5.0
```

Environment used on Niagara in `$SCRATCH/.julia-niagara`
```
  [a9b6321e] Atomix v0.1.0 ⚲
  [6a3955dd] ImageFiltering v0.7.9
  [033835bb] JLD2 v0.5.10
  [9e8cae18] Oceananigans v0.93.3 ⚲
  [d0ccf422] Oceanostics v0.14.5
  [6fe1bfb0] OffsetArrays v1.15.0
  [276daf66] SpecialFunctions v2.5.0
```

## Publication
To recreate simulations for publication, submit jobs as follows

### Mist (simulations)
`scripts/preinit/m2.sh` -> `scripts/Ri/Ri00.sh`\
`scripts/preinit/m2-Ri01.sh` -> `scripts/Ri/Ri01.sh`\
`scripts/preinit/m2-Ri02.sh` -> `scripts/Ri/Ri02.sh`\
`scripts/preinit/m4.sh` -> `scripts/amplitude/m4.sh`\
`scripts/preinit/m6.sh` -> `scripts/amplitude/m6.sh`\
`scripts/preinit/m8.sh` -> `scripts/amplitude/m8.sh`

### Niagara (post-process)
`scripts/post-process/pp-Ri00.sh`\
`scripts/post-process/pp-Ri01.sh`\
`scripts/post-process/pp-Ri02.sh`\
`scripts/post-process/pp-m4.sh`\
`scripts/post-process/pp-m6.sh`\
`scripts/post-process/pp-m8.sh`

### Figures
`fig/paper_figures.jl` creates all figures and videos included in the paper or as supplementary material. Functions that create figures are individually documented in their respective source files, located in `fig`.

## Usage
If you would like to modify the simulations, there are a few options. Running a simulation is done with\
`julia [JULIA OPTIONS] -- path/to/simulation.jl output_folder stop_time save_interval Ro Ri α β Nx Ny Nz Q ampl init_time init_rate preinit preinit_path`\
A description of each of the arguments:
### Basic arguments
 - `output_folder` Path to folder for simulation output
 - `stop_time` Simulation stop time in units of 1/f
 - `save_interval` Time between saved snapshots
 - `Ro` Minimum Rossby number of reference state
 - `Ri` Minimum Richardson number of reference state
 - `α` Ratio between filament separation and frontal jet width
 - `β` Aspect ratio, mixed layer depth relative to filament separation
 - `Nx` Number of across-front grid cells
 - `Ny` Number of down-front grid cells (size of domain will be such that horizontal grid is isotropic)
 - `Nz` Number of vertical grid cells (3/4 of which are in the mixed layer)
### Additional arguments
The remaining arguments have additional behaviour:
 - `Q` Surface cooling (flux of buoyancy upwards)
   - If zero, then the model uses a `SmagorinskyLilly` closure and `CenteredSecondOrder` advection.
   - If non-zero, then the model uses no closure and implicit `WENO` advection of order 9.
 - `ampl` Level of noise
   - If `preinit_path` is empty, then `ampl` is the maximum absolute value of the symmetric uniform distribution that initialises the `u` and `w` fields.
   - If `preinit_path` is set, then `ampl` is a prefactor that multiplies the fluctuation fields from the preinitialisation.
 - `init_time` The amount of time to restore the reference buoyancy profile for
 - `init_rate` The rate of buoyancy relaxation to the reference state during the initialisation
 - `preinit` Whether this simulation is a pre-initialisation
   - If true, the simulation runs without a filament.
   - If false, the simulation runs as normal.
 - `preinit_path` Path to pre-initialisation
   - If empty, the initial condition is some uniformly distributed noise with maximum absolute value `ampl` applied to the `u` and `w` fields.
   - If set, the initial condition is the fluctuations in the final saved timestep in `preinit_path`, multiplied by `ampl`. 
