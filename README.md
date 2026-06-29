#1 download genomes

2, use miniprot to find the rbcS 

okay I am getting hung up on distinguishing real vs false subunits. I've taken two tacks here: 
1. take ref peptides, scan genomes, get candidates 
2. get genome annotations, cluster via alignment, assign ref peptide via homology 

both necessitate some referencing for now, which is OK since I have directed hypotheses 

1. is without the biases of genome annotation, and better picks up on pseudogenes. But it is vulnerable to assembly error, and produces more results
2. is likelier to detect expressed genes, and has some more defensibility. But it is vulnerable to annotation bias, and may miss pseudogenes 

I don't think I care about pseudogenes. I think the issue with the pre-alignment in the annotation step is that it may make giant and innaccurate gene families. 

I think we want an exclusion step (not, never was, this protein), as well as a collapse step - these things are basically the same 

- exclusion - doesn't align to reference, doesn't fold to structure. 
- collapse - within X AA, with X embedding distance? 

Then when we have a clean set of genes I think we can run the association 

so the pipeline should be: 
- scan (miniprot, grabs peptides too)
- align vs ref (percent identity filter) - drop via pid/coverage (at least 50% in both)
- fold (esmfold, tm?) - drop via tm 

Keep all the sequences and embed them in plantcad and esm and do pca 
Expect to see the things the filter drops cluster distinctly 

I've set up the scan, and i've set up the alignment and fold (on rbcS, I can start here)

What I need is a tsv with 
species gene_id nucleotide_seq amino_acid_seq pid_vs_ref coverage_self coverage_ref tm_score length aln_filter_passed fold_filter_passed

Then I'll embed the sequences in esm and plantcad2 and see what is goody. 







awk -F'\t' 'NR>1 {print ">" $1 "\n" $4}' rbcS_results/rbcS_peptides.tsv | \
    /programs/mafft/bin/mafft --auto --thread 4 - > rbcS_results/rbcS_peptides_aligned.faa


    code to get species median bio8 

Rscript code/species_bio8_from_batch.R \
  --species "Quercus rubra" \
  --gbif-file /workdir/hdd29/chloroplastTempAssociation/data/0004558-260409193756587.csv \
  --bio8-tif /workdir/hdd29/chloroplastTempAssociation/data/tifs/bio08_Mean_Temperature_Wettest_Quarter.tif \
  --out data/species_bio8_medians.csv


Rscript code/species_bio8_grep.R \                                                                                                                                2 ↵
  --species-file data/species.txt \
  --gbif-file /workdir/hdd29/chloroplastTempAssociation/data/0004558-260409193756587.csv \
  --bio8-tif /workdir/hdd29/chloroplastTempAssociation/data/tifs/wc2.1_2.5m_bio_8.tif \
  --grep-engine parallel \
  --threads 10 \
  --block 500M

  Rscript code/species_bio8_grep_v3.R \
  --species-file data/genecad_species.txt \
  --gbif-file /workdir/hdd29/chloroplastTempAssociation/data/0004558-260409193756587.csv \
  --bio8-tif /workdir/hdd29/chloroplastTempAssociation/data/tifs/wc2.1_2.5m_bio_8.tif \
  --out data/genecad_species_bio8_medians.csv \
  --grep-engine parallel \
  --threads 10 \
  --block 500M \
  --force