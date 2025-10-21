#!/bin/bash
#SBATCH -n 12                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-08:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p test           # Partition to submit to
#SBATCH --mem=70G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o log/repeatmasker%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e log/repeatmasker%j.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --account=oeb275r

BFHE_assembly="cyan_flye_NP_assembly_incl_addedZandW_v2_rm1kb.fasta.gz"
RepMod_library="cyan_flye_NP_incl_addZW_v2_rm1kb-families.fa"
BFHE_assembly_RepMask="RepMask/cyan_flye_NP_assembly_incl_addedZandW_v2_rm1kb.fasta.masked"

module load RepeatMasker/4.0.5-fasrc05

# First softmask with the RepeatModeler library 
# mkdir RepMask
# RepeatMasker -xsmall -pa 12 -e ncbi -lib $RepMod_library -dir RepMask $BFHE_assembly 

# Next, use that masked assembly as input to mask with RepBase Aves library 
# mkdir RepBase_softmask 
# RepeatMasker -xsmall -pa 12 -e ncbi -lib repbase_extract_Aves.fasta -dir RepBase_softmask $BFHE_assembly_RepMask

# Mask all at once, together
# mkdir All_softmask_align
RepeatMasker -xsmall -a -pa 12 -e ncbi -lib RepMod_RepBase_cleaned_BOP_CF_libraries.fasta -dir Cleaned_all_softmask_align $BFHE_assembly
