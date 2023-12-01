#!/bin/bash
#SBATCH --job-name=AF2_predict
#SBATCH --output=AF2_predict.out
#SBATCH --error=AF2_predict.err
#SBATCH --mail-type=ALL
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20gb
#SBATCH --partition=gpu
#SBATCH --gpus=a100:1


pwd; hostname; date

module load cuda/11 gcc
export PATH="/py_path/dev/miniconda3/envs/colabfold/bin/:$PATH"
export LD_LIBRARY_PATH="/my_path/dev/miniconda3/envs/colabfold/bin/:$LD_LIBRARY_PATH"


#colabfold_batch sequence.fasta mdm2 --model-type 'alphafold2_ptm'
colabfold_batch barstar.fasta barstar --model-type 'alphafold2_ptm'
#colabfold_batch barstar_barnase.fasta barstar_barnase --model-type 'alphafold2_multimer_v3'
