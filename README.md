## Scripts associated with the manuscript "Phylogenetic distribution and recombination suppression of neo-sex chromosomes in honeyeaters (Aves: Meliphagidae)"

This readme provides an overview of the analyses conducted in this manuscript. 

The software and programs used in this manuscript include:

- [FindZX](https://github.com/hsigeman/findZX)
- [TOGA](https://github.com/hillerlab/TOGA)
  - lastz
  - NextFlow (v19.12.0)
- [flye](https://github.com/mikolmogorov/Flye)
- [NextPolish](https://github.com/Nextomics/NextPolish)
- [YAHS](https://github.com/c-zhou/yahs)
- [freebayes](https://github.com/freebayes/freebayes)
- [gblocks](https://www.biologiaevolutiva.org/jcastresana/Gblocks.html)
- [clipkit](https://github.com/JLSteenwyk/ClipKIT)
- [mafft](https://mafft.cbrc.jp/alignment/server/index.html)
- [IQTree2](https://github.com/iqtree/iqtree2)
- [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler/tree/master)
- [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker)
- [R (v4.3.2)]()
- [fastANI](https://github.com/ParBLiSS/FastANI/tree/master)

- and custom phasing scripts adapted from Sigeman, H., Ponnikas, S., Videvall, E., Zhang, H., Chauhan, P., Naurin, S., & Hansson, B. (2018). Insights into Avian Incomplete Dosage Compensation: Sex-Biased Gene Expression Coevolves with Sex Chromosome Degeneration in the Common Whitethroat. Genes, 9(8). https://doi.org/10.3390/genes9080373.


# Contents
- [Assembling genomes, polishing genomes, and scaffolding with HiC](#assembling-genomes-polishing-genomes-and-scaffolding-with-hic)
- [TOGA genome annotation](#toga-genome-annotation)
- [Phasing gametologs](#phasing-gametologs)
- [Creating alignments of loci across the neo-sex chromosomes](#creating-alignments-of-loci-across-the-neo-sex-chromosomes)
- [Expected likelihood weights](#expected-likelihood-weights)
- [Likelihood ratio test of mtDNA and nuclear topologies](#likelihood-ratio-test-of-mtDNA-and-nuclear-topologies)

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
To scaffold the genomes to chromosome-level I used the HiC assembler YAHS. First I adapted Arima pipeline to map HiC reads to the draft genomes.

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

To annotate genomes I used TOGA, which leverages whole-genome alignment with a high-quality genome (from Gallus gallus in this case) to find and annotate orthologous genes, infer exon-intron structure, gene loss, gene duplication, and identify pseudogenes.

First I created a whole-genome alignment with a [helper script](https://github.com/hillerlab/make_lastz_chains) that runs lastz with the program NextFlow. This was resource intensive and required a lot of optimization to run on a slurm-based HPC. I particularly had to optimize the size of the chunks to align for each genome: `seq1_chunk` and `seq2_chunk`. If the chunks were too large, the jobs would time-out or require too much memory. You may also have to adjust the total number of jobs that NextFlow is allowed to run and manage on the HPC.  

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
Then, I ran the TOGA python script, supplying the resulting chain file (pair-wise alignment) and the path to the query, reference, and reference annotation bed file.
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



## Phasing gametologs
## Creating alignments of loci across the neo-sex chromosomes
## Expected likelihood weights
## Likelihood ratio test of mtDNA and nuclear topologies







