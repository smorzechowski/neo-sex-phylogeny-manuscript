#!/bin/bash
#SBATCH -n 16                # Number of cores
#SBATCH -t 7-00:00           # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p bigmem            # Partition to submit to
#SBATCH --mem=250G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -J flye              # job name  
#SBATCH -o log/flye_cyan_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e log/flye_cyan_%j.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --mail-type=ALL      # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=smorzechowski@g.harvard.edu  # Email to which notifications will be sent
#SBATCH --account=oeb275r

module load python

#I had set up a Conda environment that I installed flye in. Change name as needed
source activate flye

#Replace with path to basecalled reads in FASTQ format
chry_reads="/n/holylfs04/LABS/edwards_lab/Lab/smorzechowski/meliphagid/Nanopore_data/20210908-SOrzechowski/SOrzechowski-Chry-Str/20210908_1608_2D_PAG50972_fde3b899/20210908_1608_2D_PAG50972_fde3b899_barcode11.pass.fastq.gz"
striat_reads="/n/holylfs04/LABS/edwards_lab/Lab/smorzechowski/meliphagid/Nanopore_data/20210908-SOrzechowski/SOrzechowski-Chry-Str/20210908_1608_2D_PAG50972_fde3b899/20210908_1608_2D_PAG50972_fde3b899_barcode12.pass.fastq.gz"
mel_reads="/n/holylfs04/LABS/edwards_lab/Lab/smorzechowski/meliphagid/Nanopore_data/20210907_SOrzechowski/SOrzechowski_Mel-178/20210907_1702_1B_PAG51008_7828342d/20210907_1702_1B_PAG51008_7828342d.pass.fastq.gz"
cyan_reads_sans_addedZ_NRC="/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2021-07-19/01-minimaps2/fastq_files/cyan_prom_reads_sans_addedZ_NRC.fastq.gz"
cyan_reads_sans_addedZ_NRC_updated="/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2021-07-19/01-minimaps2/fastq_files/cyan_prom_reads_sans_addedZ_NRC_updated.fastq.gz"
cyan_reads_sans_addedZ_NRC_ZPAR="/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2021-07-19/01-minimaps2/fastq_files/cyan_prom_reads_sans_addedZ_NRC_ZPAR.fastq.gz"

#--nano-raw: indicates uncorrected Nanopore reads. Options for corrected reads.
#-g: Approximate genome size.
#-o: path where assembly will be output
#--threads: number of threads used for assembly
flye --nano-raw $cyan_reads_sans_addedZ_NRC_ZPAR -g 1g -o /n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2021-09-29/01-flye/cyan_flye_assembly_sans_addedZ_NRC_ZPAR --threads 16
