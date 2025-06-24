#! /bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p shared
#SBATCH -e log/arima_hic_mapping%A.err
#SBATCH -o log/arima_hic_mapping%A.out
#SBATCH -J arima_hic
#SBATCH --mem=60g
#SBATCH --time=24:00:00
#SBATCH --account=oeb275r

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


# How to Accommodate Technical Replicates

# This pipeline is currently built for processing a single sample with a
# read1 and read2 fastq file.
# Technical replicates (eg. one library split across multiple lanes) should # be merged before running the MarkDuplicates command.
# If this step is run, the names and locations of input files to subsequent # steps will need to be modified in order for subsequent steps to run
# correctly.
# The code below is an example of how to merge technical replicates.


# REP_NUM=X # number of the technical replicate set e.g. 1
# REP_LABEL=$LABEL\_rep$REP_NUM
# INPUTS_TECH_REPS=('bash' 'array' 'of' 'bams' 'from' 'replicates') # BAM  files you want combined as technical replicates
# example bash array - INPUTS_TECH_REPS=('INPUT=A.L1.bam' 'INPUT=A.L2.bam'  'INPUT=A.L3.bam')
# java -Xmx8G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles
#$INPUTS_TECH_REPS OUTPUT=$TMP_DIR/$REP_LABEL.bam USE_THREADING=TRUE
#ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT

echo "### Step 4: Mark duplicates"
java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ \
-jar $PICARD MarkDuplicates INPUT=$PAIR_DIR/$SRA.bam OUTPUT=$REP_DIR/$REP_LABEL.bam \
METRICS_FILE=$REP_DIR/metrics.$REP_LABEL.txt TMP_DIR=$TMP_DIR \
ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

$SAMTOOLS index $REP_DIR/$REP_LABEL.bam

perl $STATS $REP_DIR/$REP_LABEL.bam > $REP_DIR/$REP_LABEL.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"
