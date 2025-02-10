#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=24:00:00
#SBATCH --job-name=SM18
#SBATCH --output=../scratch/logs/SM18_1.txt

module load julia/1.10.4
export JULIA_DEPOT_PATH=$SCRATCH/.julia-mist
export JULIA_SCRATCH_TRACK_ACCESS=0

cd ~/filament-instability-redux

STOP_TIME="18.0"
SAVE_INTERVAL="0.1"
ROSSBY="0.8"
RICHARDSON="0.001"
ALPHA="1.5"
BETA="0.1"
NX=1024
NY=128
NZ=128
COOLING="0.3"
NOISE="2e-8"
INIT_TIME="6"
INIT_RATE="1"
PREINIT="false"
PREINIT_PATH=""

julia -t 32 -- src/simulation.jl ../scratch/filament-instability-redux/cooling/SM18 $STOP_TIME $SAVE_INTERVAL $ROSSBY $RICHARDSON \
$ALPHA $BETA $NX $NY $NZ $COOLING $NOISE $INIT_TIME $INIT_RATE $PREINIT $PREINIT_PATH
