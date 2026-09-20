#!/bin/bash
#SBATCH --job-name=november_2023-sampling_completeness
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --output=november_2023-sampling_completeness-%j.out

eval "$(conda shell.bash hook)"
conda activate R-4.5

export TAR_PROJECT="november_2023-sampling_completeness"

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

echo "November 2023 sampling completeness completed successfully."
