#!/bin/bash
#SBATCH --job-name=stability_flor_pi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G                 # EDIT: 3-D FLOR fields are large; raise if the job is killed for memory
#SBATCH --time=02:00:00           # EDIT: untested on the real data, so start generous
#SBATCH --error=logs/stability_flor_pi_%j.err
#SBATCH --output=logs/stability_flor_pi_%j.out
#SBATCH --account=s1001

# Run the buoyancy-index (stability) analysis for FLOR ctrl and hitopo PI.
# Submit from the repo root:  sbatch sbatch_stability_flor.sh
# (slurm does not create the logs/ directory; the job fails silently if it is missing)

source /usr/share/modules/init/bash

# batch shells don't read the `conda init` block in ~/.bashrc, so load conda's shell functions explicitly
source /usr/local/other/anaconda/GEOSpyD/23.5.2-0_py3.11/2024-01-30/etc/profile.d/conda.sh
conda activate nam_topo

# run from the directory the job was submitted from (the repo root)
cd "${SLURM_SUBMIT_DIR}"

echo "Job ${SLURM_JOB_ID} started on $(hostname) at $(date)"
echo "Python: $(which python)"

# the script changes into notebooks/ itself, so figures land in ../figs
python notebooks/stability_flor_pi_2xco2.py

echo "Finished at $(date)"
