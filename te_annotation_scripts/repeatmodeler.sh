#!/bin/bash 
#SBATCH -n 12                # Number of cores 
#SBATCH -N 1                # Ensure that all cores are on one machine 
#SBATCH -t 7-00:00          # Runtime in D-HH:MM, minimum of 10 minutes 
#SBATCH -p shared           # Partition to submit to 
#SBATCH --mem=75G           # Memory pool for all cores (see also --mem-per-cpu) 
#SBATCH -o log/repeatmodeler%j.out  # File to which STDOUT will be written, %j inserts jobid 
#SBATCH -e log/repeatmodeler%j.err  # File to which STDERR will be written, %j inserts jobid 
#SBATCH --account=oeb275r 

MEL='/n/holylfs04/LABS/edwards_lab/Lab/smorzechowski/meliphagid'
assembly='/n/holylfs04/LABS/edwards_lab/Lab/smorzechowski/meliphagid/ReferenceAssemblies/cyan_flye_NP_assembly_incl_addedZandW_v2_rm1kb.fasta'
RepeatModeler='/n/home09/smorzechowski/bin/RepeatModeler-2.0.2a' 

# $RepeatModeler/BuildDatabase -name cyan_flye_NP_incl_addZW_v2_rm1kb -engine ncbi $assembly 
$RepeatModeler/RepeatModeler -pa 12 -engine ncbi -database cyan_flye_NP_incl_addZW_v2_rm1kb 2>&1 | tee repeatmodeler.log
