#!/bin/bash
#SBATCH -J data_analysis_for_paper_harvardrc
#SBATCH -n 32 # Number of cores requested
#SBATCH -N 1 # Number of CPUs; if 1, then requested cores will be on the same machine
#SBATCH -t 500 # Runtime in minutes
#SBATCH -p murphy # Partition to submit to; alternative is murphy_secure
#SBATCH --mem=64000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o output/data_analysis_%A_%a.out # Standard out goes to this file; %A is job name, %a is array_id for batch jobs
#SBATCH -e output/data_analysis_%A_%a.err # Standard err goes to this filehostname
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=qiantianchen@fas.harvard.edu # Edit for your own email!
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02 R/3.5.1-fasrc03 # load R
export R_LIBS_USER=$HOME/apps/R_3.5.1:$R_LIBS_USER
R CMD BATCH /n/home00/tqian/data_analysis_morogoro/data_analysis_for_paper_harvardrc.R /n/home00/tqian/data_analysis_morogoro/Rout/data_analysis_for_paper_harvardrc.Rout.${SLURM_ARRAY_TASK_ID}
