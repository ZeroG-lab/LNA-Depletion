# Description:
# This R script identifies the most abundant sequences from a set of BAM files,
# counts their occurrences, and outputs the results in a table format. 
# Additionally, it generates a heat map to visualize the relative abundance of these sequences 
# across all samples.
#
# Instructions:
# 1. Select the annotation file for your organism (preferably GTF, GFF works as well).
# 2. Set the working directory to the folder containing the BAM files.
# 3. Execute the script to obtain the output table and heat map. Run it line by line and check the output.
#
# Inputs:
# - BAM files located in the specified working directory. Must contain only one alignment per read! (for STAR use --outSAMmultNmax 1)
#
# Outputs:
# - A CSV file containing the counts of the most abundant sequences.
# - A heat map showing the abundance of the top sequences across samples.
# - A figure of merit showing expected depletion performance.
#
# Requirements:
# - R version 4.0 or higher.
#
# Authors: Dario Ricciardi and Maik Böhmer


# INSTALL AND LOAD LIBRARIES ###########################################################################

#install.packages(c("BiocManager", "dplyr","tidyr","ggplot2","viridis", "rstudioapi"))
#BiocManager::install(c("GenomicFeatures", "Rsamtools", "txdbmaker"))

# Load necessary libraries
library(GenomicFeatures)
library(txdbmaker)
library(Rsamtools)
library(rstudioapi)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)


# DEFINE FUNCTIONS AND IMPORT ANNOTATION ################################################################

# Function to count sequences and map to genes in a BAM file
count_sequences <- function(bam_file, gene_ranges) {
  
  print("Extracting sequences")
  
  # Open the BAM file
  bam <- scanBam(BamFile(bam_file))
  
  # Extract sequences and mapping positions
  seqs <- data.frame("sequences" = as.character(bam[[1]]$seq),
                     "chrom" = bam[[1]]$rname,
                     "start" = bam[[1]]$pos,
                     "end" = bam[[1]]$pos + width(DNAStringSet(as.character(bam[[1]]$seq))) - 1)
  
  print("Counting sequences and fetching gene IDs")
  
  # Count all sequences and only display uniques with their counts
  seqs.count <- count(group_by_all(seqs))
  
  # Create a GRanges object for the sequences
  sequences_gr <- GRanges(seqnames = seqs.count$chrom, ranges = IRanges(start = seqs.count$start, end = seqs.count$end), strand = "*")
  
  # Map sequences to genes
  hits <- findOverlaps(sequences_gr, gene_ranges, ignore.strand = TRUE, select = "first")
  gene_ids <- as.character(mcols(gene_ranges)[hits, "gene_id"])
  
  # Add gene ids to counted sequences
  seqs.count$gene_id <- gene_ids
  seqs.count <- data.frame(seqs.count)
  
  
  # Return data
  return(seqs.count)
}

# Function to summarize the same sequences even if they are mapped to different loci
accumulateFragments <- function(x, accumulate, ntop = 10000, fuzzyness = 3){
  contams <- data.frame()
  y <- na.omit(unique(x$Sequence[1:1000000])[1:ntop])
  y <- y[order(nchar(y))]
  ylen <- length(y)
  # Loop through all sequences from short to long
  while (length(y) > 0) {
    k <- y[1]
    #print(k) #For debugging
    # Pick out set of sequences that match the shortest one of interest (fuzzy matching if desired)
    if(accumulate == "fuzzy"){
      temp <- x[grepl(substr(k, fuzzyness+1, nchar(k)-fuzzyness), x$Sequence),]
    }else if(accumulate == "shortest"){
      temp <- x[grepl(k, x$Sequence),]
    }else if (accumulate == "sequence"){
      temp <- x[x$Sequence == k,]
    }else{print("accumulate must be either 'sequence', 'shortest' or 'fuzzy'")}
    # Sum up counts and collapse all mappings into a single character string
    temp <- temp[order(nchar(temp$Sequence)),]
    if(length(bam_files) == 1){
      temp[1,3:(2+length(sequence_counts_list))] <- sum(temp[,3:(2+length(sequence_counts_list))])
    }else{
      temp[1,3:(2+length(sequence_counts_list))] <- colSums(temp[,3:(2+length(sequence_counts_list))])
      }
    if(length(paste(unique(na.omit(temp$Gene_ID)), sep = ",")) > 0){
      temp$Gene_ID[1] <- paste(sort(unique(na.omit(temp$Gene_ID))), collapse = ",")
    }else temp$Gene_ID[1] <- NA
    
    #Bind temp to contams
    contams <- na.omit(rbind(contams, temp[1,]))
    #Eliminate matches from fragment list
    x <- x[!x$Sequence %in% unique(temp$Sequence),]
    y <- y[grep(k, y, invert = TRUE)]
    print(paste(length(y), "of", ylen, "sequences remaining"))
  }
  if(length(bam_files) == 1){
  return(contams[order(-contams[, 3:(2+length(sequence_counts_list))]),])
  }else{
    return(contams[order(-rowSums(contams[, 3:(2+length(sequence_counts_list))])),])
  }
}

# Select annotation file (GTF)
annotation_file <- file.choose() #window may appear behind R-Studio window, check if you don't see it.
# annotation_file <- "" #or define path manually

# Create a TxDb object from the annotation file using txdbmaker
txdb <- txdbmaker::makeTxDbFromGFF(annotation_file)

# Extract gene ranges
gene_ranges <- genes(txdb)


# PROCESS INPUT FILES ###################################################################################

# Set working directory containing the .bam files
setwd(selectDirectory()) #window may appear behind R-Studio window, check if you don't see it.
#setwd("") #or set manually

# Get list of BAM files in the current working directory
bam_files <- list.files(pattern = "\\.out\\.bam$")

# Initialize an empty list to store sequence counts from each file
sequence_counts_list <- list()

# Loop over each BAM file and get sequence counts
for (bam_file in bam_files) {
  print(paste("Working on", bam_file))
  sequence_counts <- count_sequences(bam_file, gene_ranges)
  samplename <- gsub("|\\.?Aligned.*bam$","", bam_file)
  sequence_counts$sample <- samplename
  sequence_counts <- sequence_counts[order(sequence_counts$n, decreasing = TRUE),]
  sequence_counts_list[[samplename]] <- sequence_counts[,c(1,5:6)]
  rm(bam_file, sequence_counts, samplename)
  print("Done")
}

# Calculate sum of all counts per sample
total_counts <- unlist(lapply(sequence_counts_list, function(x){sum(x[,2])}))


#### MERGE DATA FRAMES ###################################################################################

# Build comparison table over all samples
# Merge all data frames into one
combined_sequence_counts <- sequence_counts_list
for (i in names(combined_sequence_counts)) {
  combined_sequence_counts[[i]] <- combined_sequence_counts[[i]][,1:3]
  combined_sequence_counts[[i]]$rnd_ID <- runif(n= length(combined_sequence_counts[[i]]$gene_id))
  rm(i)
}

# Merge all lists into one (takes some time). Or in case of single file reorder columns
if(length(bam_files) == 1){
  combined_sequence_counts <- combined_sequence_counts[[1]][,c(1,3,4,2)]
  }else{
    combined_sequence_counts <- Reduce(function(x, y) merge(x, y, by = c("sequences", "gene_id", "rnd_ID"), all = TRUE), combined_sequence_counts)
    }

# Clean up random gene ids and colnames
combined_sequence_counts <- combined_sequence_counts[,-3]
colnames(combined_sequence_counts) <- c("Sequence", "Gene_ID", names(sequence_counts_list))

# Convert absolute counts to percentage of counts in relation to total counts per sample
for (i in 3:length(colnames(combined_sequence_counts))) {
 combined_sequence_counts[,i] <- combined_sequence_counts[,i] / total_counts[i-2] * 100
 rm(i)
}

# Replace NA with 0
combined_sequence_counts[is.na(combined_sequence_counts)] <- 0
combined_sequence_counts$Gene_ID[combined_sequence_counts$Gene_ID == 0] <- NA

# Sort by row-wise sums
if(length(bam_files) == 1){
  print("Skipping. Table already sorted.")
  }else{
    combined_sequence_counts <- combined_sequence_counts[order(rowSums(combined_sequence_counts[,c(names(sequence_counts_list))]), decreasing = TRUE),]
    }

# SUMMARIZE SEQUENCES ###############################################################################

#Summarize sequences across loci
# Accumulate counts of the same sequences
# "sequence" to summarize by exact sequence
# "shortest" to summarize all sequences that share a shortest common sequence
# "fuzzy" is the same as shortest but allows for overhangs of 3 nt at each end (or set fuzzyness with argument)
# Lower the ntop parameter to speed up the analysis at the cost of less optimized target sequences.
accumulated_sequences <- accumulateFragments(combined_sequence_counts, accumulate = "shortest", ntop = 10000)


# RESULTS AND PLOTTING ###############################################################################

# Write the comparison table to a CSV file
write.csv(accumulated_sequences, "./Contaminants.csv", row.names = FALSE)

# Import table for plotting
accumulated_sequences <- read.csv("./Contaminants.csv") #uncomment if you already ran the analysis and just want to plot the data or use the example dataset

# Define how many of the top contaminants you want to plot
top <- 30

# Transform data to long format for plotting and use top contaminant sequences
heatmap_data <- accumulated_sequences[1:top,] %>% #Edit columns to include all or just specific samples (e.g. [1:top, 1:4] for first 4 samples)
  gather(key = "Sample", value = "Count", -Sequence, -Gene_ID)
heatmap_data$Gene_ID <- factor(accumulated_sequences$Gene_ID[1:top], levels = unique(accumulated_sequences$Gene_ID[1:top]), ordered = TRUE)

# Add sequence length information
heatmap_data$Sequence <- paste0("...", heatmap_data$Sequence, "..."
                                #" (", nchar(heatmap_data$Sequence), ")" #uncomment to include sequence length information
                                )

# Order heatmap by average percentage
heatmap_data$Sequence <- factor(heatmap_data$Sequence, levels = rev(unique(heatmap_data$Sequence)), ordered = TRUE)

# Order heatmap by specific sample. Uncomment and insert number of sample you want to order by.
# orderSample <- unique(heatmap_data$Sample)[1]
# heatmap_data$Sequence <- factor(subset(heatmap_data, Sample == orderSample)[,1],
#                                 levels = subset(heatmap_data, Sample == orderSample)[order(subset(heatmap_data, Sample == orderSample)$Count, decreasing = FALSE),1], ordered = TRUE
#                                  )

# Calculate percentages to add to plot
if(length(bam_files) == 1){
  heatmap_data$perc[heatmap_data$Sequence == heatmap_data$Sequence[1]] <- paste0(round(sum(accumulated_sequences[1:top, 3:length(colnames(accumulated_sequences))])), "%")
}else{
  heatmap_data$perc[heatmap_data$Sequence == heatmap_data$Sequence[1]] <- paste0(round(colSums(accumulated_sequences[1:top, 3:length(colnames(accumulated_sequences))])), "%")
}

# Optional: Cleanup sample names for plotting
#heatmap_data$Sample <- gsub("PATTERN","", heatmap_data$Sample) #change PATTERN to what you want to remove from the sample name
#heatmap_data$Sample <- factor(heatmap_data$Sample, levels = unique(heatmap_data$Sample), ordered = TRUE) #uncomment to order samples in given order instead of alphabetically
#heatmap_data$Sample <- factor(heatmap_data$Sample, levels = c("SAMPLE_1", "SAMPLE_2", "ETC"), ordered =  TRUE) #uncomment and add sample names for manual ordering

# Plot contaminant heatmap
ggplot(heatmap_data, aes(x = Sample, y = Sequence, fill = Count)) +
  geom_tile() +
  scale_fill_viridis(name="Percent", option = "mako", direction = -1, begin = 0.2) +
  annotate("text",
           x = heatmap_data$Sample,
           y = top + 1,
           label = heatmap_data$perc,
           size = 2.5) +
  theme_minimal(base_size = 10)+
  theme(strip.text.y.right = element_text(angle = 0),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(),
        legend.position = "right") +
  coord_cartesian(clip = "off", ylim = c(1, top + 0.5))

# Save Plot
ggsave("./Contaminant_Heatmap.pdf",
       width = length(unique(heatmap_data$Sample))+10, height = top*0.6, units = "cm"
)


# ELBOW CURVE ANALYSIS FOR LNA SET DETERMINATION AND EXPECTED PERFORMANCE ############################

# Calculate cumulative percentages along the top contaminant fragments
heatmap_data$cumPerc <- NA
heatmap_data$topN <- NA
for (i in unique(heatmap_data$Sample)) {
  for (k in 1:nrow(subset(heatmap_data, Sample == i))) {
    heatmap_data[which(heatmap_data$Sample == i),]$cumPerc[k] <- sum(heatmap_data[which(heatmap_data$Sample == i),]$Count[1:k])
    heatmap_data[which(heatmap_data$Sample == i),]$topN[k] <- k
  }
  rm(i,k)
}

# Plot elbow curve for all samples (Rerun from "RESULTS AND PLOTTING" with higher "top" variable to see more LNA targets)
ggplot(heatmap_data, aes(x = topN, y = cumPerc, color = Sample))+
  #geom_line(linewidth = 1)+ #uncomment to see detailed graph for every sequence instead of binning
  geom_line(data = subset(heatmap_data, topN %in% c(1,seq(5,top,5))),linewidth = 1.5)+ #Binning every 5 fragments
  scale_color_viridis_d(name="Sample", option = "mako", begin = 0.1, end = 0.9) +
  labs(title = "Elbow curve of cumulative fragment percentages", x = "Number of LNAs", y = "Cumulative Contaminant Percentage")+
  scale_x_continuous(breaks = seq(1,top,1), minor_breaks = NULL)+
  scale_y_continuous(breaks = seq(0,100,10))+
  geom_vline(xintercept = seq(5,top,5))+
  theme_minimal()

# Save Plot
ggsave("./ElbowCurve.pdf",
       width = top+2, height = 16, units = "cm")
