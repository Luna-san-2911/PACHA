#!/bin/bash
#SBATCH --job-name=PACHA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=pacha-%j.log

# Generic one-allocation wrapper. Set partition, account, QoS, CPU,
# memory and time on the sbatch command line or edit the directives.
# Activate the pacha Conda environment and put PACHA on PATH before
# submission (or add your own environment activation here).
set -euo pipefail

if ! command -v PACHA >/dev/null 2>&1; then
    echo "ERROR: Activate the pacha environment and add the PACHA repository to PATH before submitting." >&2
    exit 1
fi

# All input paths and scientific/resource arguments come from the user.
# Keep local last so this wrapper never submits additional SLURM jobs.
exec PACHA "$@" --executor local
