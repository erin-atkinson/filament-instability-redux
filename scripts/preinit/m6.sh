#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=1:00:00
#SBATCH --job-name=preinit-m6
#SBATCH --output=../scratch/logs/filament-instability.preinit.m6.txt

module load julia/1.10.4
export JULIA_DEPOT_PATH=$SCRATCH/.julia-mist
export JULIA_SCRATCH_TRACK_ACCESS=0

cd ~/filament-instability-redux

STOP_TIME="1.0"
SAVE_INTERVAL="1.0"
ROSSBY="0.8"
RICHARDSON="0.01"
ALPHA="1.5"
BETA="0.1"
NX=1024
NY=1024
NZ=128
COOLING="0"
NOISE="2e-6"
INIT_TIME="0"
INIT_RATE="0"
PREINIT="true"
PREINIT_PATH=""

julia -t 32 -- src/simulation.jl ../scratch/filament-instability-redux/preinit/m6 $STOP_TIME $SAVE_INTERVAL $ROSSBY $RICHARDSON \
$ALPHA $BETA $NX $NY $NZ $COOLING $NOISE $INIT_TIME $INIT_RATE $PREINIT $PREINIT_PATH
