#! /bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p test
#SBATCH -e log/gphase_I%A.err
#SBATCH -o log/gphase_I%A.out
#SBATCH -J gphase_I
#SBATCH --mem=5g
#SBATCH --time=12:00:00


####################################
# Author: Hanna Sigeman, 2018
# Contact: hanna.sigeman@biol.lu.se
#
# This script is the first of two used to phase whole-genome sequence data into a Z and W chromosome gene sequences using one female and one male sample.
# The script uses a VCF file where genotypes from one female and one male have been called for each genomic position within exon regions of specified genes.
#
# Usage: ./genotype_phasing_script_I.sh <VCF file> <output>
#
# Comment: <VCF file> ($1) should specify a VCF file where the first sample (in column 10) is a male sample (having a ZZ genotype) and the
# second sample (in column 11) is a female (having a ZW genotype). <output> ($2) is the chosen output name for the modified VCF file
# (which will be used as input in the second script).
#
# Output: A modified VCF file where each position is marked with the allele corresponding to the Z and W chromosome.
#
# Modified: Sophie Orzechowski, January 2024
# Modifications:
# Only flagging depth (DP) up to 10
# Move depth line to the top of if else statements, since homomorphic sites may have a low depth due to deletions


cat $1 | grep -v "^#" | sed -E "s/([ATCG]+),([ATCG]+)/\1\t\2/1" | awk --re-interval '{
if($6 ~ /^[A-Z]+,[A-Z]+/) print $0,"ZW","multiN",length($4)
else if($6 ~ /^[A-Z]+/ && $5!=$6) print $0
else print $1,$2,$3,$4,$5,$5,$6,$7,$8,$9,$10,$11}' | tr ' ' '\t' | awk -v OFS="\t" --re-interval '{
if (/multi/) print $0
else if(/DP=[0-9];/) print $0,"ZW","depth10",length($4)
else if ($5==".") print $0,"ZW","same",length($4)
else if ($11=="." || $12==".") print $0,"ZW","missing",length($4)
else if ($7<=20) print $0,"ZW","qual20",length($4)
else if ($11 ~ /^0\/0:/ && $12 ~ /^0\/0:/) print $0,"Z","ref",length($4),"break",$0,"W","ref",length($4)
else if ($11 ~ /^0\/1:/ && $12 ~ /^0\/0:/) print $0,"Z","N",length($4),"break",$0,"W","N",length($4)
else if ($11 ~ /^0\/0:/ && $12 ~ /^0\/1:/) print $0,"Z","ref",length($4),"break",$0,"W","alt1",length($4)
else if ($11 ~ /^0\/1:/ && $12 ~ /^0\/1:/) print $0,"Z","N",length($4),"break",$0,"W","N",length($4)
else if ($11 ~ /^1\/1:/ && $12 ~ /^0\/1:/) print $0,"Z","alt1",length($4),"break",$0,"W","ref",length($4)
else if ($11 ~ /^0\/1:/ && $12 ~ /^1\/1:/) print $0,"Z","N",length($4),"break",$0,"W","N",length($4)
else if ($11 ~ /^1\/1:/ && $12 ~ /^1\/1:/) print $0,"Z","alt1",length($4),"break",$0,"W","alt1",length($4)
else if ($11 ~ /^1\/2:/ && $12 ~ /^1\/1:/) print $0,"Z","N",length($4),"break",$0,"W","alt1",length($4)
else if ($11 ~ /^1\/1:/ && $12 ~ /^1\/2:/) print $0,"Z","alt1",length($4),"break",$0,"W","alt2",length($4)
else if ($11 ~ /^1\/2:/ && $12 ~ /^1\/2:/) print $0,"Z","N",length($4),"break",$0,"W","N",length($4)
else if ($11 ~ /^2\/2:/ && $12 ~ /^1\/2:/) print $0,"Z","alt2",length($4),"break",$0,"W","alt1",length($4)
else if ($11 ~ /^1\/2:/ && $12 ~ /^2\/2:/) print $0,"Z","N",length($4),"break",$0,"W","N",length($4)
else if ($11 ~ /^2\/2:/ && $12 ~ /^2\/2:/) print $0,"Z","alt2",length($4),"break",$0,"W","alt2",length($4)
else print $0,"ZW","Un",length($4)}' | tr ' ' '\t' | sed 's/break/\n/g' | sed -e 's/^[ \t]*//' | sed -e 's/[ \t]*$//' | grep -v ">" | while read a b c d e f g h i j k l m n o ; do printf $a"\t"$b"\t"$c"\t"$d"\t"$e"\t"$f"\t"$g"\t"$h"\t"$i"\t"$j"\t"$k"\t"$l"\t"$m"\t"$n"\t"$o"\t" ; eval printf -- 'N%.s' {1..$o} ; echo ; done > $2
