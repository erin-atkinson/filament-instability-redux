#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=12:00:00
#SBATCH --job-name=Ri00
#SBATCH --output=../scratch/logs/filament-instability.Ri00.txt

module load julia/1.10.4
export JULIA_DEPOT_PATH=$SCRATCH/.julia-mist
export JULIA_SCRATCH_TRACK_ACCESS=0

cd ~/filament-instability-redux

STOP_TIME="40.0"
SAVE_INTERVAL="0.1"
ROSSBY="0.8"
RICHARDSON="0.01"
ALPHA="1.5"
BETA="0.1"
NX=1024
NY=1024
NZ=128
COOLING="0"
NOISE="1"
INIT_TIME="0"
INIT_RATE="0"
PREINIT="false"
PREINIT_PATH="../scratch/filament-instability-redux/preinit/m2"

julia -t 32 -- src/simulation.jl ../scratch/filament-instability-redux/Ri/Ri00 $STOP_TIME $SAVE_INTERVAL $ROSSBY $RICHARDSON \
$ALPHA $BETA $NX $NY $NZ $COOLING $NOISE $INIT_TIME $INIT_RATE $PREINIT $PREINIT_PATH
