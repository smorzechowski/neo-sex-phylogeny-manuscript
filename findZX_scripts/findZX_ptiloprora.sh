#/bin/bash

# Example of how to run FindZX snakemake pipeline
# Create a file of chromosomes from bTaeGut1.4.pri to highlight in the findZX plots!

nano chromosomes_highlight.txt
NC_044241.2
NC_044217.2
NC_044212.2

# LOGIN node:  holylogin04
tmux new -s new_sesh

conda activate findingZX
snakemake -s workflow/findZX-synteny -j 8 -R all --configfile Ptiloprora/config.yml --cluster-config Ptiloprora/cluster.yml \
--cluster " sbatch -p shared -t 07-00:00 -n 12 --mem 40G -N 1 -A oeb275r "  --use-conda

# COMPLETED successfully! Now create the html doc!

tmux new -s ptiloprora_report

#started an srun session from holylogin04
conda activate findingZX

snakemake -s workflow/findZX-synteny --configfile Ptiloprora/config.yml --cores 1 -R all -k --use-conda --report report_Ptiloprora.html
