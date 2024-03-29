# RAS16S
Repository of pilot microbial study

### MiSeq analysis:

**Step 1:**
For the initial processing of the fasta files: **DADA2_Combo.R**
Takes as an input the demultiplexed fasta files. Outputs excel tables with ASVs and taxa

**Step2:**
Post-processing of the data obtained from Step 1 (Phyloseq.xlsx): **Pyloseq_Combo.R**
- Summarizes the contents of the data
- Downsamples Reads to Even Depth
- Removes Singletons
- Computes the total and average prevalances of the features of each phylum
- Removes phyla with only 1 Hit and NAs
- Calculates Alpha diversity
- Calculates Alpha observed (richness)
- Calculates Alpha Shannon diversity
- Tests for significant difference in species diversity
- Plots the Alpha diversities
- Plots MDS
- Performs two-way PERMANOVA

### PacBio analysis: 

**dada2_PacBio_GitHub.R** 
