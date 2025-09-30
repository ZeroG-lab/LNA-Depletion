# Repository for the publication titled: Data-driven design of LNA-blockers for efficient contaminant removal in Ribo-seq libraries
## Abstract
Ribo-Seq libraries often contain abundant non-coding RNA contaminants, which, because of their high sequence variability and diverse fragmentation, are challenging to remove. We present a computational pipeline that identifies experiment-specific target sequences and allows for their efficient depletion using custom LNA probes in a single pipetting step, thereby increasing sequencing yield and reducing costs. A public LNA repository will support sharing validated targets within the research community.

<img width="500" height="500" alt="GraphicalAbstract" src="https://github.com/user-attachments/assets/88297401-e86a-4c1d-9361-ee79d4f1908c" />

## Contaminant Profiles
Summarized and ordered contaminant tables for *Arabidopsis thaliana* across all tested growth conditions can be found in /ContaminantTables.

## Instructions
To run the script on your own data, please first make sure that you have aligned your reads in such a way that there is EXACTLY one alignment per read.
For example, when using STAR provide the flag --outSAMmultNmax 1.

Acquire an annotation file for your organism (GTF or GFF).

Open the Contam_Ident.R script in R-Studio with a current R version (4.0+).

Install and load all necessary packages as indicated in the script and follow the included instructions line by line.

In case you just want to test the script, skip ahead to line 202 and import the provided Contaminants.csv file.
This file contains the summarized contaminants from our growth condition tests and lets you try out the plotting functions.

**Attention**: The script currently expects at least two samples. If you only have a single sample, make a copy of the .bam file and pretend to have two samples. A single sample will lead to errors. This problem will be fixed in the future.
