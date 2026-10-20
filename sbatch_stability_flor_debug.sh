#!/bin/bash
#SBATCH --job-name=stability_flor_debug
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G                 # EDIT: smaller than the full run since only 10 years are loaded
#SBATCH --time=00:30:00           # EDIT: keep within the debug QoS time limit
#SBATCH --qos=debug
#SBATCH --error=logs/stability_flor_debug_%j.err
#SBATCH --output=logs/stability_flor_debug_%j.out
#SBATCH --account=s1001

# Quick DEBUG version of sbatch_stability_flor.sh: short wall time, debug QoS, and only 10 PI years.
# Submit from the repo root:  mkdir -p logs && sbatch sbatch_stability_flor_debug.sh
# Estimate the wait without submitting:  sbatch --test-only sbatch_stability_flor_debug.sh
#
# --- Even faster: debug interactively (no requeue between attempts) ---
# 1) get a shell on a compute node (from the login node):
#      salloc --ntasks=1 --cpus-per-task=2 --mem=32G --time=00:30:00 --qos=debug --account=s1001
# 2) once the prompt returns on the compute node, from the repo root:
#      source /usr/local/other/anaconda/GEOSpyD/23.5.2-0_py3.11/2024-01-30/etc/profile.d/conda.sh
#      conda activate nam_topo
#      export STAB_YR_MIN=0101 STAB_YR_MAX=0110 STAB_FIG_DIR=../figs/debug
#      python notebooks/stability_flor_pi_2xco2.py        # rerun after each fix; no queue wait
#      exit                                               # release the allocation when done

source /usr/share/modules/init/bash

# batch shells don't read the `conda init` block in ~/.bashrc, so load conda's shell functions explicitly
source /usr/local/other/anaconda/GEOSpyD/23.5.2-0_py3.11/2024-01-30/etc/profile.d/conda.sh
conda activate nam_topo

# run from the directory the job was submitted from (the repo root)
cd "${SLURM_SUBMIT_DIR}"

# read by notebooks/stability_flor_pi_2xco2.py: only 10 years, figures to a separate debug folder
export STAB_YR_MIN=0101
export STAB_YR_MAX=0110
export STAB_FIG_DIR=../figs/debug     # relative to notebooks/ (the script changes into it)

echo "Job ${SLURM_JOB_ID} started on $(hostname) at $(date)"
echo "Python: $(which python)"
echo "Years: ${STAB_YR_MIN}-${STAB_YR_MAX}; figures -> ${STAB_FIG_DIR}"

python notebooks/stability_flor_pi_2xco2.py

echo "Finished at $(date)"
