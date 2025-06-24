#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p edwards,shared
#SBATCH -e log/ssphase_II%A.err
#SBATCH -o log/gphase_II%A.out
#SBATCH -J ssphase_II
#SBATCH --mem=5g
#SBATCH --time=2-00:00

####################################
# Author: Hanna Sigeman, 2018
# Contact: hanna.sigeman@biol.lu.se
#
# This script is the second of two used to phase whole-genome sequence data into a single haplotype using one male sample.
# The script uses a modified VCF file which is the output of script number 1. It also requires a 5 column bed file (0-based positions)
# where the first three columns specifies the genome ranges of the exons, the fourth column contains the name of the gene, and the fifth the cds number
#
# Usage: ./single_sample_cds_phasing_script_II.sh <BED file> <VCF file> <output>
#
# Comment:
# <BED file> ($1) is a bed file containing genomic ranges for exons/cds. HAS to be pre-sorted by cds (from 0 to x)
# <VCF file> ($2) is the output from script 1. <BED file> ($1) is a bed file containing genomic ranges for exons/cds.
# <output> ($3) is the chosen output name for the fasta file containing the phased gene sequences.
# Output: A fasta file containing haplogype where the sequence headers correspond to the gene names specified in the bed file ($1).
#
# Modified: Sophie Orzechowski, January 2024
#
# removed  awk '{if($1=="'"$cont"'") print $0}' from third line; condition already satisfied
# added +0 to awk statement to force numerical format of intervals; otherwise there's a bug
# changed depth20 to depth10
# removed || $14=="covDiff" and $14=="sizediffN" from last line
# adjusted the column numbers to reflect the # columns in the first phasing script -- hopefully correctly!
#
#####################################

species=$4

cat $1  | awk -v OFS="\t" '{print $1,$2,$3,$4,$5}' | while read cont start end gene cds
do echo ">${species}_${gene}_${cont}_cds${cds}|${start}:${end}" ; \
cat $2 | grep -v "^#" | awk '{if($2>="'"$start"'"+0 && $2<="'"$end"'"+0) print $0}' | awk '{

# Select haplotypes for single sample (male)
if ($12=="M") print $0}' | awk --re-interval '{
if ($13=="ref" || $13=="same") print $4
else if ($13=="alt") print $5
else if ($13=="alt1") print $5
else if ($13=="alt2") print $6
else if ($13=="N" || $13=="missing" || $13=="qual20" || $13=="depth10" || $13=="multiN" || $13=="Un") print $15}' | tr ' ' '\t' | tr -d "\n"  ; done | sed 's/>/\n>/' > $3

