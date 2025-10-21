#!/bin/bash
#SBATCH -n 6                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 7-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared           # Partition to submit to
#SBATCH --mem=75G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o log/nextpolish%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e log/nextpolish%j.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --account=oeb275r

#1. Prepare sgs_fofn
# ls reads_R1.fq reads1_R2.fq > sgs.fofn

#2. Create run.cfg
# genome=input.genome.fa
# echo -e "task = best\ngenome = $genome\nsgs_fofn = sgs.fofn" > run.cfg


#3. Run
/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2021-09-29/04-nextpolish/NextPolish/nextPolish ./cyanotis/run.cfg
/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2021-12-07/01-nextpolish/NextPolish/nextPolish ./mel/run.cfg
