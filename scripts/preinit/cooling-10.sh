#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=4:00:00
#SBATCH --job-name=preinit-strong-cooling-10
#SBATCH --output=../scratch/logs/filament-instability.preinit.strong-cooling-10.txt

module load julia/1.10.4
export JULIA_DEPOT_PATH=$SCRATCH/.julia-mist
export JULIA_SCRATCH_TRACK_ACCESS=0

cd ~/filament-instability-redux

STOP_TIME="10.0"
SAVE_INTERVAL="1.0"
ROSSBY="0.8"
RICHARDSON="0.001"
ALPHA="1.5"
BETA="0.1"
NX=1024
NY=1024
NZ=128
COOLING="5e-2"
NOISE="1e-8"
INIT_TIME="0"
INIT_RATE="0"
PREINIT="true"
PREINIT_PATH=""

julia -t 32 -- src/simulation.jl ../scratch/filament-instability-redux/preinit/stronger-cooling-10 $STOP_TIME $SAVE_INTERVAL $ROSSBY $RICHARDSON \
$ALPHA $BETA $NX $NY $NZ $COOLING $NOISE $INIT_TIME $INIT_RATE $PREINIT $PREINIT_PATH
