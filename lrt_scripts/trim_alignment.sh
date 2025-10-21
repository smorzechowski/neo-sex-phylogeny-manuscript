#!/bin/bash

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

