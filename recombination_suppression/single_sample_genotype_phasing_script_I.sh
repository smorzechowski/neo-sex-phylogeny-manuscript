#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p test
#SBATCH -e log/ss_phase_I%A.err
#SBATCH -o log/ss_phase_I%A.out
#SBATCH -J ssphase_I
#SBATCH --mem=5g
#SBATCH --time=12:00:00


cat $1 | grep -v "^#" | sed -E "s/([ATCG]+),([ATCG]+)/\1\t\2/1" | awk --re-interval '{
if($6 ~ /^[A-Z]+,[A-Z]+/) print $0,"M","multiN",length($4)
else if($6 ~ /^[A-Z]+/ && $5!=$6) print $0
else print $1,$2,$3,$4,$5,$5,$6,$7,$8,$9,$10}' | tr ' ' '\t' | awk -v OFS="\t" --re-interval '{
if (/multi/) print $0
else if(/DP=[0-9];/) print $0,"M","depth10",length($4)
else if ($5==".") print $0,"M","same",length($4)
else if ($11==".") print $0,"M","missing",length($4)
else if ($7<=20) print $0,"M","qual20",length($4)
else if ($11 ~ /^0\/0:/) print $0,"M","ref",length($4)
else if ($11 ~ /^0\/1:/) print $0,"M","N",length($4)
else if ($11 ~ /^1\/1:/) print $0,"M","alt1",length($4)
else if ($11 ~ /^1\/2:/) print $0,"M","N",length($4)
else if ($11 ~ /^2\/2:/) print $0,"M","alt2",length($4)
else print $0,"M","Un",length($4)}' | tr ' ' '\t' | sed -e 's/^[ \t]*//' | sed -e 's/[ \t]*$//' | grep -v ">" | \
while read a b c d e f g h i j k l m n ; do printf $a"\t"$b"\t"$c"\t"$d"\t"$e"\t"$f"\t"$g"\t"$h"\t"$i"\t"$j"\t"$k"\t"$l"\t"$m"\t"$n"\t" ; \
eval printf -- 'N%.s' {1..$n} ; echo ; done > $2
############################################################################################################
