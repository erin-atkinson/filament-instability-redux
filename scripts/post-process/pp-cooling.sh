#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=4:00:00
#SBATCH --job-name=pp-strong-cooling
#SBATCH --output=../scratch/logs/filament-instability.post-process.strong-cooling.txt

module load julia/1.10.4
export JULIA_DEPOT_PATH=~/.julia-niagara
export JULIA_SCRATCH_TRACK_ACCESS=0
# Path to ramdisk to avoid too many writes on parallel filesystem
mkdir /dev/shm/filament-instability-pp
export RAM=/dev/shm/filament-instability-pp

cd ~/filament-instability-redux

# Location of output.jld2
export SIM_OUTPUT_FOLDER=../scratch/filament-instability-redux/cooling/strong-cooling

julia -t 40 -- src/analysis/post-process.jl $SIM_OUTPUT_FOLDER DFM $RAM
julia -t 40 -- src/analysis/post-process.jl $SIM_OUTPUT_FOLDER PV $RAM
julia -t 40 -- src/analysis/post-process.jl $SIM_OUTPUT_FOLDER TKE $RAM
julia -t 40 -- src/analysis/post-process.jl $SIM_OUTPUT_FOLDER PSI $RAM
julia -t 40 -- src/analysis/post-process.jl $SIM_OUTPUT_FOLDER FRONT $RAM
