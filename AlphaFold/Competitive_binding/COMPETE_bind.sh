#!/bin/bash
#SBATCH --job-name=compete_bind
#SBATCH --output=compete.out
#SBATCH --error=compete.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ymartineznoa@ufl.edu
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20gb
#SBATCH --account=alberto.perezant
#SBATCH --partition=gpu
#SBATCH --gpus=a100:1


pwd; hostname; date

module load cuda/11 gcc
export PATH="/blue/alberto.perezant/liweichang/dev/miniconda3/envs/colabfold/bin/:$PATH"
export LD_LIBRARY_PATH="/blue/alberto.perezant/liweichang/dev/miniconda3/envs/colabfold/bin/:$LD_LIBRARY_PATH"


colabfold_batch sequences.fasta MDMX_MDM2 --model-type 'alphafold2_ptm'
