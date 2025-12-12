#!/bin/bash
#SBATCH --job-name=duckweed_and_box
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --output=duckweed_and_box-%j.out

eval "$(conda shell.bash hook)"
conda activate R-4.5

export TAR_PROJECT="duckweed_and_box"

# execute the tar_make command
if ! Rscript -e "targets::tar_make()"; then
    echo "Error: tar_make command failed."
    exit 1
fi

# prune the targets
if ! Rscript -e "targets::tar_prune()"; then
    echo "Error: tar_prune command failed."
    exit 1
fi
