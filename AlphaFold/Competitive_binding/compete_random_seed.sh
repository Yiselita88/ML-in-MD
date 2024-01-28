# #!/bin/bash
#SBATCH --job-name=compete_bind
#SBATCH --output=compete.out
#SBATCH --error=compete.err
#SBATCH --mail-type=ALL
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20gb
#SBATCH --partition=gpu
#SBATCH --gpus=a100:1


pwd; hostname; date

module load cuda/11 gcc
export PATH="/your_path/dev/miniconda3/envs/colabfold/bin/:$PATH"
export LD_LIBRARY_PATH="/your_path/dev/miniconda3/envs/colabfold/bin/:$LD_LIBRARY_PATH"


colabfold_batch sequences.fasta result_num-ensemble --model-type 'alphafold2_ptm' --num-seeds 1 --random-seed 1