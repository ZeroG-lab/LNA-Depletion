# Repository for the publication titled: Data-driven design of LNA-blockers for efficient contaminant removal in Ribo-Seq libraries
Published in Scientific Reports [https://doi.org/10.1038/s41598-026-43117-3](https://doi.org/10.1038/s41598-026-43117-3)
## Abstract
Ribo-Seq libraries often contain highly abundant non-coding RNA contaminants, which are challenging to remove due to their high sequence variability and diverse fragmentation patterns. We present an organism-independent computational pipeline that identifies experiment-specific target sequences and enables their efficient depletion using custom-tailored LNA probes in a single pipetting step. We demonstrate that LNA-based depletion is most effective during library amplification and has no effect on gene-level quantification. Contaminant depletion in Arabidopsis libraries nearly doubled the yield of coding reads, significantly improving cost-effectiveness.

<img width="500" height="500" alt="GraphicalAbstract" src="https://github.com/user-attachments/assets/88297401-e86a-4c1d-9361-ee79d4f1908c" />

## Contaminant Profiles
Summarized and ordered contaminant tables for *Arabidopsis thaliana* across all tested growth conditions can be found in /ContaminantTables.

## Instructions
To run the script on your own data, please first make sure that you have aligned your reads in such a way that there is EXACTLY one alignment per read.
For example, when using STAR provide the flag --outSAMmultNmax 1.

Acquire an annotation file for your organism (GTF or GFF).

Open the Contam_Ident.R script in R-Studio with a current R version (4.0+).

Install and load all necessary packages as indicated in the script and follow the included instructions line by line.

If you just want to test the script, use the alignments and annotation supplied in the ExampleData folder (unzip annotation before using).
The alignments consist of two *Arabidopsis thaliana* Ribo-Seq libraries. One of them undepleted (No_LNA_1) and the other depleted with the five LNAs listed in the manuscript (LNA_1).
