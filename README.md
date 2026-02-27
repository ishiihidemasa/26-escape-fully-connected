This repository contains Julia codes and simulation datasets.  

# Reproducing figures from simulation results
To reproduce figures based on simulation results stored in `data` directory, 
the easiest way is:  

- download this repository
- open the directory in VSCode
- 'Reopen in Container' using `Dev Container`
  - requires Docker Desktop
- Done! You are ready to run (most) Julia scripts in this repository!

Alternatively, you can create the appropriate Julia environment locally -- not in Docker container --
using `Project.toml` and `Manifest.toml` in the root directory.

# Reproducing simulation results
We performed direct numerical simulations of model SDEs on a high-performance computer (HPC) using `apptainer`.  
Roughly speaking, the procedure was as follows:

- Prepare an apptainer container on a HPC where you can run `run-simulation-global.jl`.
  - We used the following shell script to initialize a container according to `Project.toml` and `Manifest.toml`.
 ```
#!/bin/sh
apptainer build julia.sif docker://julia:1.11.4-bookworm
apptainer exec julia.sif julia --project=. -e "using Pkg; Pkg.instantiate()"
```
- Run `run-simulation-global.jl` in the container, passing parameter values of your choise.
  - You might also want to modify parameters defined in the Julia script.
  - I note that multiprocessing is used: be aware of CPU and memory consumptions!

Please contact me (ISHII Hidemasa) if you need further clarification.
