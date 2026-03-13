#!/bin/bash

#SBATCH --job-name=rgh_simulation
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=2000
#SBATCH --gres=disk:5000
#SBATCH -c 1
#SBATCH --time=24:00:00
#SBATCH --mail-user=%u@jlab.org

$RGH_SIM_HOME/job.sh
