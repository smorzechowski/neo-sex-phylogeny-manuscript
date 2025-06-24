## Scripts associated with the manuscript "Phylogenetic distribution and recombination suppression of neo-sex chromosomes in honeyeaters (Aves: Meliphagidae)

This readme provides an overview of the analyses conducted in this manuscript. 

The software and programs used in this manuscript include:

- [FindZX](https://github.com/hsigeman/findZX)
- [TOGA](https://github.com/hillerlab/TOGA)
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
- [Assembling genomes and scaffolding with HiC](#assembling-genomes-and-scaffolding-with-hic)
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
To 


## TOGA genome annotation
## Phasing gametologs
## Creating alignments of loci across the neo-sex chromosomes
## Expected likelihood weights
## Likelihood ratio test of mtDNA and nuclear topologies







