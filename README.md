# filament-instability-redux
Setup for a

Simulations and post processing done on SciNet Mist, Niagara systems respectively.

## Publication
To recreate figures for publication, submit jobs as follows

### Mist (simulations)
`scripts/preinit/m2.sh` -> `scripts/Ri/Ri00.sh`\
`scripts/preinit/m2-Ri01.sh` -> `scripts/Ri/Ri01.sh`\
`scripts/preinit/m2-Ri02.sh` -> `scripts/Ri/Ri02.sh`\
`scripts/preinit/m4.sh` -> `scripts/amplitude/m4.sh`\
`scripts/preinit/m6.sh` -> `scripts/amplitude/m6.sh`\
`scripts/preinit/m8.sh` -> `scripts/amplitude/m8.sh`\
`scripts/preinit/cooling-10.sh`

### Niagara (post-process)
`scripts/post-process/pp-Ri00.sh`\
`scripts/post-process/pp-Ri01.sh`\
`scripts/post-process/pp-Ri02.sh`\
`scripts/post-process/pp-m4.sh`\
`scripts/post-process/pp-m6.sh`\
`scripts/post-process/pp-m8.sh`\
`scripts/post-process/pp-cooling-preinit.sh`

### Figures
`paper_figures.jl` creates all figures\


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
