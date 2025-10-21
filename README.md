# Scripts associated with the manuscript "Neo-sex chromosomes track the mitochondrial phylogeny and exhibit an extensive added stratum of recombination suppression in honeyeaters (Aves: Meliphagidae)"

This readme provides an overview of the analyses conducted in this manuscript. 

The software and programs used in this manuscript include:

- [FindZX](https://github.com/hsigeman/findZX)
- [TOGA](https://github.com/hillerlab/TOGA)
  - lastz
  - NextFlow (v19.12.0)
- [flye](https://github.com/mikolmogorov/Flye)
- [NextPolish](https://github.com/Nextomics/NextPolish)
- [YAHS](https://github.com/c-zhou/yahs)
- [BWA (v.7.17)](https://bio-bwa.sourceforge.net/bwa.shtml)
    - picard 
    - samtools (v1.20)
- [freebayes](https://github.com/freebayes/freebayes)
- [mafft](https://mafft.cbrc.jp/alignment/server/index.html)
- [gblocks](https://www.biologiaevolutiva.org/jcastresana/Gblocks.html)
- [clipkit](https://github.com/JLSteenwyk/ClipKIT)
- [IQTree2](https://github.com/iqtree/iqtree2)
- [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler/tree/master)
- [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker)
- [R (v4.3.2)]()
- [fastANI](https://github.com/ParBLiSS/FastANI/tree/master)

- and custom phasing scripts adapted from Sigeman, H., Ponnikas, S., Videvall, E., Zhang, H., Chauhan, P., Naurin, S., & Hansson, B. (2018). Insights into Avian Incomplete Dosage Compensation: Sex-Biased Gene Expression Coevolves with Sex Chromosome Degeneration in the Common Whitethroat. Genes, 9(8). https://doi.org/10.3390/genes9080373.


# Contents
- [Assembling genomes, polishing genomes, and scaffolding with HiC](#assembling-genomes-polishing-genomes-and-scaffolding-with-hic)
- [TOGA genome annotation](#toga-genome-annotation)
- [Calling variants and phasing gametologs](#calling-variants-and-phasing-gametologs)
- [Creating alignments of loci across the neo-sex chromosomes](#creating-alignments-of-loci-across-the-neo-sex-chromosomes)
- [Expected likelihood weights](#expected-likelihood-weights)
- [Likelihood ratio test of mtDNA and nuclear topologies](#likelihood-ratio-test-of-mtDNA-and-nuclear-topologies)
- [Scanning for neo-sex chromosomes with FindZX](#scanning-for-neo-sex-chromosomes-with-findzx)

## Assembling genomes, polishing genomes, and scaffolding with HiC

To assemble the genomes with Oxford Nanopore long reads, I used the assembler flye v.2.8.1:
```
#--nano-raw: indicates uncorrected Nanopore reads. Options for corrected reads.
#-g: Approximate genome size.
#-o: path where assembly will be output
#--threads: number of threads used for assembly

flye --nano-raw $cyan_reads_sans_addedZ_NRC_ZPAR -g 1g \
-o /n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2021-09-29/01-flye/cyan_flye_assembly_sans_addedZ_NRC_ZPAR --threads 16
```

To polish genomes, I used the program NextPolish: 
```
cat run.cfg

[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 2
multithread_jobs = 5
genome = ./mel_flye_assembly.fasta
genome_size = auto
workdir = ./01_rundir_trimmed
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = ./sgs.fofn
sgs_options = -max_depth 100

[lgs_option]
lgs_fofn = ./lgs.fofn
lgs_options = -min_read_len 5k -max_depth 100
lgs_minimap2_options = -x map-ont


#Run
/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2021-12-07/01-nextpolish/NextPolish/nextPolish ./mel/run.cfg
```
To scaffold the genomes to chromosome-level I used the HiC assembler YAHS. First I adapted an Arima pipeline to map HiC reads to the draft genomes.

```
source activate samtools
module load jdk/20.0.1-fasrc01

SRA='cyan-204-874802_S3HiC_combined'
LABEL='HiC_cyan_204'
BWA='/n/home09/smorzechowski/bin/bwa/bwa'
SAMTOOLS='samtools'
#IN_DIR='/n/holylfs04/LABS/edwards_lab/Lab/smorzechowski/meliphagid/HiC_data'
IN_DIR='/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2023-06-14/01-yahs/data'
REF='/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2023-06-14/01-yahs/cyan_flye_NP_assembly_incl_addedZandW_v2_rm1kb.fasta'
FAIDX='$REF.fai'
PREFIX='bwa_index'
RAW_DIR='/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2023-06-14/01-yahs/raw'
FILT_DIR='/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2023-06-14/01-yahs/filtered'
FILTER='/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2023-06-14/01-yahs/filter_five_end.pl'
COMBINER='/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2023-06-14/01-yahs/two_read_bam_combiner.pl'
STATS='/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2023-06-14/01-yahs/get_stats.pl'
PICARD='/n/home09/smorzechowski/bin/picard/build/libs/picard.jar'
TMP_DIR='/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2023-06-14/01-yahs/tmp'
PAIR_DIR='/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2023-06-14/01-yahs/paired/'
REP_DIR='/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2023-06-14/01-yahs/dedup/'
REP_LABEL=$LABEL\_rep1
MERGE_DIR='/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2023-06-14/01-yahs/replicates'
MAPQ_FILTER=10
CPU=12

echo "### Step 0: Check output directories exist & create them as needed"
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
[ -d $REP_DIR ] || mkdir -p $REP_DIR
[ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR


# Run only once! Skip this step if you have already generated BWA index files
#echo "### Step 0: Index reference"
#$BWA index -a bwtsw $REF

#echo "### Step 1.A: FASTQ to BAM (1st)"
#$BWA mem -t $CPU $REF $IN_DIR/$SRA\_R1.fastq.gz | $SAMTOOLS view -@ $CPU -Sb  - > $RAW_DIR/$SRA\_R1.bam

#echo "### Step 1.B: FASTQ to BAM (2nd)"
#$BWA mem -t $CPU $REF $IN_DIR/$SRA\_R2.fastq.gz | $SAMTOOLS view -@ $CPU -Sb  - > $RAW_DIR/$SRA\_R2.bam

#echo "### Step 2.A: Filter 5' end (1st)"
#$SAMTOOLS view -h $RAW_DIR/$SRA\_R1.bam | perl $FILTER | $SAMTOOLS view  -Sb - > $FILT_DIR/$SRA\_R1.bam

#echo "### Step 2.B: Filter 5' end (2nd)"
#$SAMTOOLS view -h $RAW_DIR/$SRA\_R2.bam | perl $FILTER | $SAMTOOLS view  -Sb - > $FILT_DIR/$SRA\_R2.bam

#echo "### Step 3A: Pair reads & mapping quality filter"
#perl $COMBINER $FILT_DIR/$SRA\_R1.bam $FILT_DIR/$SRA\_R2.bam $SAMTOOLS $MAPQ_FILTER | $SAMTOOLS view -bS -t $FAIDX - | $SAMTOOLS sort -@ $CPU -o $TMP_DIR/$SRA.bam -

#echo "### Step 3.B: Add read group"
#java -Xmx4G -Djava.io.tmpdir=temp/ -jar $PICARD AddOrReplaceReadGroups INPUT=$TMP_DIR/$SRA.bam OUTPUT=$PAIR_DIR/$SRA.bam ID=$SRA LB=$SRA SM=$LABEL PL=ILLUMINA PU=none

echo "### Step 4: Mark duplicates"
java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ \
-jar $PICARD MarkDuplicates INPUT=$PAIR_DIR/$SRA.bam OUTPUT=$REP_DIR/$REP_LABEL.bam \
METRICS_FILE=$REP_DIR/metrics.$REP_LABEL.txt TMP_DIR=$TMP_DIR \
ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

$SAMTOOLS index $REP_DIR/$REP_LABEL.bam

perl $STATS $REP_DIR/$REP_LABEL.bam > $REP_DIR/$REP_LABEL.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"
```
Next, I ran YAHS to scaffold the contigs using the HiC-aligned reads.
```
source activate pretextmap
# installed samtools within pretextmap env too!
module load jdk/20.0.1-fasrc01


out="cyan_out"
outdir="."
contigs="cyan_flye_NP_assembly_incl_addedZandW_v2_rm1kb.fasta" # need to be indexed, i.e., ${test}.contigs.fasta.gz.fai is presented
hicaln="/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2023-06-14/01-yahs/dedup/HiC_cyan_204_rep1.bam" # could be .bed, .bam or .bin file

#### run yahs scaffolding
/n/home09/smorzechowski/bin/yahs/yahs -o ${outdir}/${out} ${contigs} ${hicaln} || exit 1

```

## TOGA genome annotation

To annotate genomes I used TOGA, which leverages whole-genome alignment with a high-quality genome (from *Gallus gallus* in this case) to find and annotate orthologous genes, infer exon-intron structure, gene loss, gene duplication, and identify pseudogenes.

First I created a whole-genome alignment with a [helper script](https://github.com/hillerlab/make_lastz_chains) that runs lastz with the program NextFlow. This was resource intensive and required a lot of optimization to run on a slurm-based HPC. I particularly had to optimize the size of the chunks to align for each soft-masked genome: `seq1_chunk` and `seq2_chunk`. If the chunks were too large, the jobs would time-out or require too much memory. You may also have to adjust the total number of jobs that NextFlow is allowed to run and manage on the HPC.  

```
#! /bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p edwards
#SBATCH -e log/chains%A.err
#SBATCH -o log/chains%A.out
#SBATCH -J chains
#SBATCH --mem=20g
#SBATCH --time=30-00:00
#SBATCH --account=edwards_lab

module load python
source activate lastz_chain_dependencies
module purge
#module load jdk/20.0.1-fasrc01

export PERL5LIB=/


export PATH=/n/home09/smorzechowski/bin/nextflow-19.12.0-edge:$PATH
export PATH=/n/home09/smorzechowski/bin/make_lastz_chains/kent_binaries:$PATH
export PATH=/n/home09/smorzechowski/bin/make_lastz_chains/GenomeAlignmentTools/src:$PATH

which python
which java

# Genomes must be soft-masked

target_name=$1
query_name=$2
target_genome=$3
query_genome=$4
project_directory=$5

python /n/home09/smorzechowski/bin/make_lastz_chains/make_chains.py ${target_name} ${query_name} ${target_genome} ${query_genome} \
--executor slurm --executor_partition edwards --executor_queuesize 200 --project_dir ${project_directory} --chaining_memory 100000 \
--cluster_parameters '-A edwards_lab --time=15-00:00' --seq1_chunk 30000000 --seq2_chunk 7000000

# --force_def
#--continue_arg 'lastz' If you need to restart at this point
#--continue_arg 'chainRun' if you need to restart at this point
# 'chainRun'
# 'cat' needed more memory for full genome, increasing from 4 GB to 40 GB, finally worked...
# 'chainRun' to fix tmp directory location
```
Then, I ran the TOGA python script, supplying the resulting chain file (pairwise alignment) and the path to the query, reference, and reference annotation bed file. I ran TOGA on the entire Blue-faced Honeyeater genome, but I extracted the contigs/chromosomes homologous to Chr Z and Chr 5 in Noisy Miner, Striated Pardalote, Helmeted Honeyeater, Yellow-rumped Thornbill. 
```
module load python

which java
export PATH=/n/home09/smorzechowski/bin/nextflow-19.12.0-edge:$PATH
export PATH=/n/home09/smorzechowski/bin/make_lastz_chains/kent_binaries:$PATH

TOGAPATH='/n/home09/smorzechowski/bin/TOGA'
CHAINFILE=$1
ANNOBED=$2
REFPATH=$3
QUERYPATH=$4
ISOFORMS=$5
NFCONFIG=$6
PROJECT=$7
#U12=$7

python $TOGAPATH/toga.py $CHAINFILE $ANNOBED ${REFPATH} ${QUERYPATH} --kt --pn $PROJECT -i $ISOFORMS --nc ${NFCONFIG} --cb 1,5,15,50,100 --cjn 500 --ms
```
## Calling variants and phasing gametologs

I used custom scripts adapted from Sigeman et al. 2018 to phase W-linked and Z-linked gametologs from loci on the neo-sex chromosomes. 

The scripts involve phasing variants from VCF files of male and female samples using rules of expected hetero- and homozygosity on sex chromosomes when mapping each sex a homogametic (ZZ) reference genome.  

To create VCF files, I made a homogametic version of the Blue-faced Honeyeater genome by removing the neo-W and then mapped Illumina short reads from a male and female Blue-faced Honeyeater, White-eared Honeyeater, Little Friarbird, which all possess neo-sex chromosomes. I used BWA, picard, and samtools to create bam files, mark duplicates, add read groups, validate bam files.  

```
# Input parameters
RAW_DIR=$1
DEDUP_DIR=$2
TMP_DIR=$3
GENOME=$4
INDEXBASE=$5
FASTQ1=$6
FASTQ2=$7
CPU=$8
FLOWCELL=$9
LANE=${10}
FLOWCELL_BARCODE=${11}
SAMPLE=${12}
PLATFORM=${13}
LIBPREP=${14}

module purge
module load python
source activate samtools
module load jdk/20.0.1-fasrc01

PICARD='/n/home09/smorzechowski/bin/picard/build/libs/picard.jar'
# this comes from ArimaGenomics mapping_pipeline
STATS='/n/home09/smorzechowski/bin/get_stats.pl'
BWA_PATH='/n/home09/smorzechowski/bin/bwa'

#Check output directories exist & create them as needed
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $DEDUP_DIR ] || mkdir -p $DEDUP_DIR

# create bwa index if necessary - only do this once
[ -f ${INDEXBASE}.sa ] || $BWA_PATH/bwa index -p $INDEXBASE $GENOME

# call bwa and convert to bam, directly sort
# https://www.biostars.org/p/319730/
# -M : "mark shorter split hits as secondary
# -R : Read group information
# '@RG\tID:$ID\tSM:$SM\tLB:$LB\tPU:$PU\tPL:$PL'

$BWA_PATH/bwa mem $INDEXBASE $FASTQ1 $FASTQ2 -t $CPU -M -R $(echo "@RG\tID:$FLOWCELL"_"$LANE\tSM:$SAMPLE\tLB:$LIBPREP\tPU:$FLOWCELL_BARCODE"_"$LANE"_"$SAMPLE\tPL:$PLATFORM") \
| samtools sort -@ $CPU --output-fmt BAM -o $RAW_DIR/${SAMPLE}_bwa_sorted.bam

# index the sorted bam files
samtools index $RAW_DIR/${SAMPLE}_bwa_sorted.bam

# get summary of the alignment
samtools flagstat $RAW_DIR/${SAMPLE}_bwa_sorted.bam > $RAW_DIR/${SAMPLE}_bwa_sorted.bam.flagstat

# Mark duplicates with Picard
#-XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ \
java -Xmx30G \
-jar $PICARD MarkDuplicates \
TMP_DIR=$TMP_DIR \
INPUT=$RAW_DIR/${SAMPLE}_bwa_sorted.bam \
OUTPUT=$DEDUP_DIR/${SAMPLE}_bwa_sorted_dedup.bam \
METRICS_FILE=$DEDUP_DIR/${SAMPLE}_dedup_metrics.txt \
REMOVE_DUPLICATES=false \
TAGGING_POLICY=All \
ASSUME_SORTED=true \

# Index the dedup bam files again
samtools index $DEDUP_DIR/${SAMPLE}_bwa_sorted_dedup.bam

# output stats on mapping using ArimaGenomics perl script
perl $STATS $DEDUP_DIR/${SAMPLE}_bwa_sorted_dedup.bam > $DEDUP_DIR/${SAMPLE}_bwa_sorted_dedup.bam.stats

# validate the bam files
java -Xmx30G \
-jar $PICARD ValidateSamFile \
      I=$DEDUP_DIR/${SAMPLE}_bwa_sorted_dedup.bam \
      MODE=SUMMARY
```

To call variants, I used the parallel version of program FreeBayes and used the flag `--report-monomorphic` to include invariant sites for future phylogenetic alignment purposes. 

```
BEDFILE=$1
GENOME=$2
BAMLIST=$3
SPECIES=$4
REGIONS=$5

module load python
source activate freebayes

export PERL5LIB=/

# call variants assuming a diploid sample and report monomorphic (all sites, even those that aren't variable)
freebayes-parallel $REGIONS 20 -f ${GENOME} --bam-list ${BAMLIST} --report-monomorphic > ${SPECIES}_neoZ_Ecyan_HiC_v1.0_sansW_cds_longest_isoform.vcf

# remove duplicates in the vcf file and save to new file
sbatch bcftools.jobscript ${SPECIES}_neoZ_Ecyan_HiC_v1.0_sansW_cds_longest_isoform
```
With the VCF files in hand for each species, I ran a custom phasing script adapted from Sigeman et al 2018. See methods for details on how each site was called as Z-linked or W-linked reference or alternate, etc by comparing variant patterns in the male and female sample. 

See genotype_phasing_script_I.sh in files. 

```
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

```

The above file created an annotated VCF file with flags for each site, that the second phasing script collected and turned into a fasta file of phased W and Z sequences for each locus based on a bed file of genomic ranges of CDS for each locus.

See Wlinked_cds_phasing_script_II.sh and Zlinked_cds_phasing_script_II.sh in files. 

```
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

```

To include other honeyeaters in my alignments where I only had a single male sample, including Painted Honeyeater *(Grantiella picta)* and White-naped Honeyeater *(Melithreptus lunatus)*, I adapted the phasing scripts to phase just Z-linked sequences across the neo-sex chromosomes.

See single_sample_genotype_phasing_script_I.sh in files.

```
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
```
Then I extracted the phased sequences for the single haplotype.

See single_sample_cds_phasing_script_II.sh.

```
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

```
## Creating alignments of loci across the neo-sex chromosomes

To include sequences from TOGA annotated genomes of Helmeted Honeyeater, Noisy Miner, Striated Pardalote, Yellow-rumped Thornbill, I extracted CDS sequences and subset the loci on Chr Z and Chr 5.

See agat_get_cds_fasta.jobscript in files. 

```
# Example: get fasta for Chry
gff=/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2023-11-28/04-annotation/galGal6_to_Chry_5Z.annotation.gtf
genome=/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2023-11-28/04-annotation/Chry_Chr5_Z_Ragtag.fasta
fold $genome > folded_Chry_Chr5_Z_Ragtag.fasta
output=Chry_Chr5_Z_transcripts.fasta
genome_fold=folded_Chry_Chr5_Z_Ragtag.fasta
sbatch agat_get_cds_fasta.jobscript $gff $genome_fold $output

gff=$1
genome=$2
output=$3

export PERL5LIB=/

singularity exec --cleanenv /n/home09/smorzechowski/bin/agat_0.8.0--pl5262hdfd78af_0.sif \
agat_sp_extract_sequences.pl --gff $gff --fasta $genome -o $output -t cds --plus_strand_only
```

This is an example of subsetting transcripts for one species.
```
# subset transcripts for Striat, headclean
sequences='/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2023-11-28/18-alignments/Striat_Chr5_Z_transcripts.fasta'
awk '{print $1}' $sequences > headcleaned_Striat_Chr5_Z_transcripts.fasta
genes='Ecyan_transcript_id_present_longest_transcript.txt'
sed -e 's/\.[[:digit:]]\+\.[[:digit:]]\+//g' $genes > genes_stripped.txt
grep -f genes_stripped.txt headcleaned_Striat_Chr5_Z_transcripts.fasta | awk '{print $1}' | sed -e 's/^>//g' > gene_list.txt
sbatch seqtk.jobscript headcleaned_Striat_Chr5_Z_transcripts.fasta gene_list.txt

```

I also added the species name to each fasta sequence.

```
sed 's/>/>Striat_/I' subset_headcleaned_Striat_Chr5_Z_transcripts.fasta > namefix_subset_headcleaned_Striat_Chr5_Z_transcripts.fasta
sed 's/>/>Chry_/I' subset_headcleaned_Chry_Chr5_Z_transcripts.fasta > namefix_subset_headcleaned_Chry_Chr5_Z_transcripts.fasta
sed 's/>/>Mmel_/I' subset_headcleaned_Mmel_Chr5_Z_transcripts.fasta > namefix_subset_headcleaned_Mmel_Chr5_Z_transcripts.fasta
sed 's/>/>HeHo_/I' subset_headcleaned_HeHo_Chr5_Z_transcripts.fasta > namefix_subset_headcleaned_HeHo_Chr5_Z_transcripts.fasta

```

Then I split the resulting fasta files into a separate file for each locus.

See split_fasta_files.jobscript in files. 

```
fasta=$1
species=$2


# First, split each subset fasta file into separate files per sequence
# Do this in a separate results directory I guess

cat $fasta | awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fasta")}
        print $0 > filename
}'
```

Next, for the species on which I had run the custom phasing scripts, I had to combine all phased exons into a single line of sequence for each locus. In actuality, these were phased CDS regions because I was phasing based on a TOGA annotation of Blue-faced Honeyeater and TOGA only aligns and retains CDS information, not entire exons.

Note that it was very important to concatenate the separate CDS sequences in the right order, otherwise the final transcript would be garbled and in the wrong order to be able to align to other species where I extracted transcripts from a genome annotation. Reverse complementing to get the opposite strand would not help in this case. It was essential that CDS 1, CDS 2, CDS 3, etc. were concatenated in the right order. 

See combine_exons.jobscript in files.

```
bed=$1
fasta=$2
haplotype=$3
species=$4
output=$5

#don't think I want to sort the bed file, just keep unique transcript names

###
awk '{print $4}' $bed | uniq | while read G; do cat $fasta | paste  - -  | grep -F "${haplotype}_${G}" \
| tr "\t" "\n" | awk -v G=$G -v hap=$haplotype -v spec=$species 'NR==1 {printf(">%s_%s_%s\n",spec,hap,G);} (NR%2==0) {print}' >> $output; done

###

# Now remove the extra line breaks:
cat $output | awk '!/^>/ { printf "%s", $0; n = "\n" }  /^>/ { print n $0; n = "" }  END { printf "%s", n } ' > formatted_${output}
```

I did the same for the species where I had a single male sample -- see combine_exons_single_sample.jobscript in files.


Next, I split the fasta files into separate files per sequence and then concatenated all files of the same transcript together (to make a multi-species fasta file for phylogenetic alignment.

See format_fasta_alignments.jobscript in files.

```
#fasta=$1
#species=$2
#transcripts=$3
#filepath=/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2023-11-28/18-alignments/results

# First, split each subset fasta file into separate files per sequence
# Do this in a separate directory I guess

#cat $fasta | awk '{
#        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fasta")}
#        print $0 > filename
#}'

# Next, combine (concatenate) all files of the same transcript together

# strip the transcript names
#sed -e 's/\.[[:digit:]]\+\.[[:digit:]]\+//g' $transcripts > transcripts_stripped.txt

#for transcript in $(< transcripts_stripped.txt);
#do cat ${species}_*${transcript}*.fasta > ${transcript}_alignment.fasta; done

```

To align the phased sequences, I used mafft.

```
module load python
source activate te_annot


# Round 1: alignment the phased sequences
#echo "the fasta file is $1"
#fileprefix=`echo $1 |sed 's/\.fasta//g'`
#mafft ${fileprefix}.fasta > ${fileprefix}.maf
```

To add the extracted sequences from the species with TOGA annotations, I used mafft with the flag `--adjustdirectionaccurately` to ensure that the correct strand was aligned to the rest of the sequences.

```
# Round 2: add additional sequences in the right order
echo "the alignment file is $1"
newseq=$2
fileprefix=`echo $1 |sed 's/\.maf//g'`
#basenameprefix=basename $1

file=${1##*/}
basenameprefix=`echo $file |sed 's/\.maf//g'`

mafft --add $newseq --genafpair --adjustdirectionaccurately --maxiterate 1000 --nuc --reorder ${fileprefix}.maf > ${basenameprefix}_added.maf
```

***MOST importantly, I manually inspected each alignment file to determine that the sequences had been added in the right order at all steps along the way!***

I used gblocks to trim the alignments.

```
#!/bin/bash

module load python
source activate gblocks

Gblocks $alignment -t=c -p=y -e=.gbc  
Gblocks $alignment -t=d -p=y -e=.gbd
```

## Expected likelihood weights

Once I had the final, trimmed alignment files for each locus across the neo-sex chromosomes, it was very straightforward to run the expected likelihood weight analyses in IQTree. I manually created tree topologies representing each scenario of recombination suppression against a nuclear tree backbone and mitochondrial tree backbone. I did this in the program Mesquite. Then, I provided IQTree with these predefined topologies and tested each alignment against them all with a suite of metrics, including expected likelihood weights.

```
module load python
source activate iqtree

alignment=$1
trees=$2
output=$3

# Run ELW analysis
iqtree -s $alignment --trees $trees --test 100000 --test-au --prefix $output

#-redo
```

I used awk and sed to combine all results together and create a csv file that I could manipulate in R. 

```

cat NewPAR_Ecyan_transcript_id_present_longest_transcript.txt | sed -e 's/\.[[:digit:]]\+\.[[:digit:]]\+//g' | sed -e 's|$|_alignment_revsort_added_renamed_all.maf.gbc.iqtree|' > NewPAR_iqtree_files.txt
cat $(cat NewPAR_iqtree_files.txt) > NewPAR.iqtree.reports 

cat NewZ_Ecyan_transcript_id_present_longest_transcript.txt | sed -e 's/\.[[:digit:]]\+\.[[:digit:]]\+//g' | sed -e 's|$|_alignment_revsort_added_renamed_all.maf.gbc.iqtree|' > NewZ_iqtree_files.txt
cat $(cat NewZ_iqtree_files.txt) > NewZ.iqtree.reports

cat OldZ_Ecyan_transcript_id_present_longest_transcript.txt | sed -e 's/\.[[:digit:]]\+\.[[:digit:]]\+//g' | sed -e 's|$|_alignment_revsort_added_renamed_all.maf.gbc.iqtree|' > OldZ_iqtree_files.txt
cat $(cat OldZ_iqtree_files.txt) > OldZ.iqtree.reports

# gather results
awk '/USER TREES/{flag=1;next}/ALISIM COMMAND/{flag=0}flag' NewZ.iqtree.reports > NewZ.iqtree.logtable.reports
awk '/USER TREES/{flag=1;next}/ALISIM COMMAND/{flag=0}flag' NewPAR.iqtree.reports > NewPAR.iqtree.logtable.reports
awk '/USER TREES/{flag=1;next}/ALISIM COMMAND/{flag=0}flag' OldZ.iqtree.reports > OldZ.iqtree.logtable.reports

```

I manually finished formatting the files and visualized results in ELW_recombination_suppression_figures.R -- see files. 

## Likelihood ratio test of mtDNA and nuclear topologies

I conducted a likelihood ratio test of the mtDNA and nuclear species tree topologies to confirm the presence of mitonuclear discordance. I used data from Andersen, M. J., McCullough, J. M., Friedman, N. R., Peterson, A. T., Moyle, R. G., Joseph, L., & Nyári, Á. S. (2019). Ultraconserved elements resolve genus-level relationships in a major Australasian bird radiation (Aves: Meliphagidae). Emu - Austral Ornithology, 119(3), 218–232.

I manually created trees in Mesquite matching the mtDNA and nuclear topologies in Andersen et al. 2019. I randomly subset 1000 loci from the supermatrix of UCE loci used to build the species tree. I used the partition file to extract a subset of the loci. 

```
shuf -n 1000 andersen_2019_UCE_partitions_sorted.txt > andersen_2019_UCE_partitions_1000_subsampled.txt

# and sort again
awk -F'=' '{ split($2, r, "-"); gsub(/[^0-9]/, "", r[1]); print r[1] "\t" $0 }' andersen_2019_UCE_partitions_1000_subsampled.txt | sort -n | cut -f2- > andersen_2019_UCE_partitions_1000_subsampled_sorted.txt

# extract site ranges and build a site list

awk -F= '{ gsub(/[^0-9\-]/, "", $2); split($2, r, "-"); for (i = r[1]; i <= r[2]; i++) print i }' andersen_2019_UCE_partitions_1000_subsampled_sorted.txt > andersen_2019_UCE_partitions_1000_subsampled_site_list.txt


# remap the partition file - important!
awk -F= 'BEGIN {pos=1}
{
    split($2, r, "-")
    len = r[2] - r[1] + 1
    printf("%s = %d-%d\n", $1, pos, pos+len-1)
    pos += len
}' andersen_2019_UCE_partitions_1000_subsampled_sorted.txt > andersen_2019_UCE_partitions_1000_subsampled_sorted_remapped.txt

# trim full alignment
awk 'NR==1 {
    num_taxa = $1
    aln_len = $2
    next
}
{
    taxa[NR-1] = $1
    seqs[NR-1] = $2
}
END {
    while ((getline line < "andersen_2019_UCE_partitions_1000_subsampled_site_list.txt") > 0) pos[line] = 1

    new_len = 0
    for (j = 1; j <= aln_len; j++) {
        if (pos[j]) new_len++
    }

    print num_taxa, new_len

    for (i = 1; i <= num_taxa; i++) {
        out = ""
        for (j = 1; j <= length(seqs[i]); j++) {
            if (pos[j]) out = out substr(seqs[i], j, 1)
        }
        print taxa[i], out
    }
}' genera59_mafft_min75percent_gblocks0.65_clean_nexus.phylip > genera59_mafft_min75percent_gblocks0.65_clean_nexus_trimmed_alignment.phylip

```

Then I ran IQTree to test the two topologies 
```
module load python
source activate iqtree

# Full dataset
alignment='genera59_mafft_min75percent_gblocks0.65_clean_nexus.phylip'
partitions='andersen_2019_UCE_partitions_sorted.txt'
trees='Andersen_2019_treeset.txt'

# Subsampled 25% dataset
#alignment='genera59_mafft_min75percent_gblocks0.65_clean_nexus_trimmed_alignment.phylip'
#partitions='andersen_2019_UCE_partitions_1000_subsampled_sorted_remapped.txt'
#trees='Andersen_2019_treeset.txt'

# match PartitionFinder - test partition merging only
#iqtree -s $alignment -p $partitions -z $trees -m TESTMERGEONLY -zb 10000 -zw -au -pre topology_test_multithread -T 8

# perform model selection, merging, branch optimization before topology testing
iqtree -s $alignment -p $partitions -z $trees -m MFP+MERGE -zb 10000 -zw -au -pre topology_test_full_dataset -T 12
```

## Scanning for neo-sex chromosomes with FindZX

I used the handy snakemake pipeline FindZX from Hanna Sigeman to scan for neo-sex chromosomes using genomic signatures of sex-linkage from whole-genome resequencing data: sex differences in heterozygosity (indicative of young sex chromosomes), and sex differences in mapping depth (indicative of old sex chromosomes). Importantly, FindZX requires a homogametic reference genome to map resequencing reads from males and females. I found that it is possible to use just a single individual of each sex, with fairly low coverage (5-10x), to reliably detect neo-sex chromosomes with these metrics. 

Below is an example of how I used FindZX. I ran the pipeline in a tmux terminal session so that I could do other things while it was running. The FindZX session in tmux would submit jobs to the cluster as each rule was successfully finished. It is possible to customize resource and time SLURM requests for each submitted job with the cluster.yml file. However, I ended up overriding this with the `--cluster` command where I made the amount of time, threads, and memory uniformly the same for all submitted jobs. I couldn't justify spending the extra time to refine the amount of resources requested for the many rules.   

Here is an example config.yml file:
```
## ================================= ##
## findZX config file (test dataset) ##
## ================================= ##

# Variables marked with "[findZX]" or "[findZX-synteny]" are only used when deploying 
# the pipeline with either snakefile. Other variables are used for all analyses. 

threads_max: 13
mem_max: 40000 
  # Specify maximum number of cores [threads_max] and memory [mem_max] allocation


# ============================ #
# Analysis name and input data #

run_name: Ptiloprora 
  # Select an analysis name. Output files will be stored under "results/[run_name]"

units: ./Ptiloprora/units.tsv
  # Path to sample information file

ref_genome: /n/holylfs04/LABS/edwards_lab/Lab/smorzechowski/meliphagid/ReferenceAssemblies/BFHE_MedakaAssembly/BFHE_Medaka_V1_rm1kb_refined_addedZ_V2.fasta   
  # Path to study-species reference genome (not .gz format)

synteny_ref: /n/holylfs04/LABS/edwards_lab/Lab/smorzechowski/meliphagid/ZebraFinchAssembly/GCF_003957565.2_bTaeGut1.4.pri/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna
  # [findZX-synteny] Path to synteny-species reference genome (not .gz format)

synteny_name: T_guttatus
  # [findZX-synteny] Synteny-species name (can be any string, will be used for file and directory names)


# ================= #
# Plotting settings #

window_sizes: [25000, 50000, 100000]
  # Choose genome window sizes for plotting (as many as you want)
  # Optimal sizes depend on reference genome fragmentation and size of the sex-linked region
  # Recommended sizes to start are: [50000, 100000, 1000000] (i.e. 50 kb, 100 kb, 1Mb)

chr_file: ./Ptiloprora/contigs.txt
  # [findZX] Specify a file with list of chromosomes to only plot these (otherwise leave as "None")

chr_highlight: ["ONT_addedZ_V2_contig_4_segment0"]
  # [findZX] Specify chromosomes/scaffolds to highlight in plot type 4, or leave empty

synteny_chr_file: ./Ptiloprora/chromosomes.txt
  # [findZX-synteny] Specify a file with list of chromosomes to only plot these (otherwise leave as "None")

synteny_chr_highlight: ["NC_044217.2"]
  # [findZX-synteny] Specify chromosomes/scaffolds to highlight in plot type 4, or leave empty


# ================================== #
# Trimming and subsampling of reads  #

## These three variables control trimming and subsampling of reads
## Set all to "false" to disable trimming and subsampling
## Only one variable is allowed to be "true"

trim_reads: true 
  # Set to true for trimming of reads

trim_and_subsample: false
  # Set to true for trimming and subsampling of reads

subsample_only: false
  # Set to true for subsampling of reads (but not trimming)

subsample_basepairs: 1888226
  # Specify the total number of basepairs to extract from both fastq files
  # Will be used if [trim_and_subsample] or [subsample_basepairs] is set to "true"

  # Use this script to calculate expected coverage:
  # ./code/subsampling_cov_calv.sh <REF.fasta> <WANTED_COV> 


# ========================== #
# findZX-specific parameters # 
# ===== (edit if needed) ===== #
 
mismatch_settings: [0.0, 0.2]
  # Genome coverage results will be generated from the original BAM files ("unfiltered"),
  # and two other (modifiable) mismatches settings.
  # "0.0" = 0 mismatches allowed
  # "0.2" = <=2 mismatches allowed

minSizeScaffold: "10000"
  # The mimimum size of scaffolds in the reference genome to be included in the results


# ============================ #
# External software parameters # 
# ===== (edit if needed) ===== #

params:
  trimmomatic:
  # Control Trimmomatic settings here
    pe:
      trimmer:
        - "LEADING:3"
        - "TRAILING:3"
        - "SLIDINGWINDOW:4:15"
        - "MINLEN:36"
        - "ILLUMINACLIP:workflow/meta/adapters/TruSeq3-PE.fa:2:30:10"
```

Here is how I run the snakemake pipeline: 

```
# Create a file of chromosomes from bTaeGut1.4.pri to highlight in the findZX plots!  
nano chromosomes_highlight.txt 
NC_044241.2  
NC_044217.2 
NC_044212.2

# LOGIN node:  holylogin04
tmux new -s new_sesh

conda activate findingZX 
snakemake -s workflow/findZX-synteny -j 8 -R all --configfile Ptiloprora/config.yml --cluster-config Ptiloprora/cluster.yml --cluster " sbatch -p shared -t 07-00:00 -n 12 --mem 40G -N 1 -A oeb275r "  --use-conda

# COMPLETED successfully! Now create the html doc! 

tmux new -s ptiloprora_report 

started an srun session from holylogin04 
conda activate findingZX 

snakemake -s workflow/findZX-synteny --configfile Ptiloprora/config.yml --cores 1 -R all -k --use-conda --report report_Ptiloprora.html
```





