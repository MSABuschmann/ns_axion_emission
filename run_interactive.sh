salloc --nodes 1 --qos interactive --time 02:00:00 --constraint cpu --account=m3166

srun shifter --image=docker:mil1/povray:latest --volume="/global/cfs/cdirs/m3166/buschman/QCD/povray:/temp" povray /temp/ani.pov
srun shifter --image=docker:mil1/povray:latest --volume="/global/cfs/cdirs/m3166/buschman/QCD/povray:/temp" /bin/bash

#!/bin/bash
#SBATCH --constraint=cpu
#SBATCH --nodes=1
##SBATCH --image=docker:mil1/povray:latest
##SBATCH --volume="/global/homes/b/buschman/cori/GPU/AMREX_Sims/ET-Integration/Exec/postevolution_Higgs_lH_povray/image:/temp"
#SBATCH --qos=interactive
#SBATCH --time=01:00:00
#SBATCH -A m3166

#srun shifter povray /temp/ani.pov
