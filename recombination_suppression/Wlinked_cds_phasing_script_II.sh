#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p edwards,shared
#SBATCH -e log/gphase_II%A.err
#SBATCH -o log/gphase_II%A.out
#SBATCH -J gphase_II
#SBATCH --mem=5g
#SBATCH --time=2-00:00


####################################
# Author: Hanna Sigeman, 2018
# Contact: hanna.sigeman@biol.lu.se
#
# This script is the second of two used to phase whole-genome sequence data into a Z and W chromosome gene sequences using one female and one male sample.
# The script uses a modified VCF file which is the output of script number 1. It also requires a 4 column bed file (0-based positions)
# where the first three columns specifies the genome ranges of the exons and the fourth column contains the name of the gene.
#
# Usage: ./genotype_phasing_script_II.sh <BED file> <VCF file> <output>
#
# Comment: <VCF file> ($2) is the output from script 1. <BED file> ($1) is a bed file containing genomic ranges for exons or cds.
# <output> ($3) is the chosen output name for the fasta file containing the phased gene sequences.
#
# Output: A fasta file containing phased Z and W gene sequences where the sequence headers correspond to the gene names specified in the bed file ($1).
#
# Modified: Sophie Orzechowski, January 2024
# bed file is a 5 column pre-sorted file of all cds from 0 to x -- need to remake this since I deleted it :(
# remove this awk '{if($1=="'"$cont"'") print $0}' from third line -- already know this condition is satisfied
# remove $14=="sizediffN" || $14=="covDiff" from last line -- those flags are not in the phased file
#####################################


species=$4

cat $1 | awk -v OFS="\t" '{print $1,$2,$3,$4,$5}' | while read cont start end gene cds
do echo ">${species}_Wlinked_${gene}_${cont}_${cds}|${start}:${end}" ; \
cat $2 | grep -v "^#" | awk '{if($2>="'"$start"'"+0 && $2<="'"$end"'"+0) print $0}' | awk '{

# Select W or common alleles
if ($13=="W" || $13=="ZW") print $0}' | awk --re-interval '{
if ($14=="ref" || $14=="same") print $4
else if ($14=="alt") print $5
else if ($14=="alt1") print $5
else if ($14=="alt2") print $6
else if ($14=="N" || $14=="missing" || $14=="qual20" || $14=="depth10" || $14=="multiN" || $14=="Un") print $16}' | tr ' ' '\t' | tr -d "\n"  ; done | sed 's/>/\n>/' > $3
