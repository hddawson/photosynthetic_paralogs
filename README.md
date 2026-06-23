#1 download genomes

2, use miniprot to find the rbcS 


awk -F'\t' 'NR>1 {print ">" $1 "\n" $4}' rbcS_results/rbcS_peptides.tsv | \
    /programs/mafft/bin/mafft --auto --thread 4 - > rbcS_results/rbcS_peptides_aligned.faa