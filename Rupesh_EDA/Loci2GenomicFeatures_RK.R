#!/usr/bin/env Rscript

###################################################################################
# File: Loci2GenomicFeatures_RK.R
# Input: bed (3 column bed file witout a header)
# output: directory with 4 outputs (2 plots png and pdf) and annotated bed & TSV
# Author: Rupesh Kesharwani
# Last update: Aug 30, 2024
# Copyright (c) 2024 Kesharwani RK
###################################################################################

##########################
##    Hackathon 2024    ##
##    BCM, Houston      ##
##        TDB           ##
##########################

## https://github.com/collaborativebioinformatics/tandemrepeats/tree/main

options(showArgsList = FALSE)
suppressWarnings(library(argparse))

# Set up argument parser
parser <- ArgumentParser(description = "Annotate BED file with genomic features and generate distribution plots")
parser$add_argument("-b", "--bed", required = TRUE, help = "Path to input BED file")
parser$add_argument("-t", "--tsv", required = TRUE, help = "Path to input TSV file")
parser$add_argument("-n", "--name", required = TRUE, help = "Filename prefix for output files")
parser$add_argument("-o", "--output", required = TRUE, help = "Directory to save the output files")

# Parse arguments
args <- parser$parse_args()

# Function to load R packages
load_packages <- function(packages) {
  invisible(lapply(packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    suppressPackageStartupMessages(suppressWarnings(library(pkg, character.only = TRUE, quietly = TRUE)))
  }))
}

# List of packages to load
packages <- c(
  "GenomicRanges",
  "rtracklayer",
  "GenomicFeatures",
  "ggplot2",
  "ggpubr",
  "GenomicDistributions",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "snpStats",
  "GenomicTools.fileHandler",
  "GenomicDistributionsData"
)

# Load the packages
load_packages(packages)

# Load the BED file
bed_file <- args$bed
bed_gr <- import(bed_file, format = "BED")

# Load a TxDb object (use your appropriate genome)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Extract genomic features
exons <- exons(txdb)
utr5 <- fiveUTRsByTranscript(txdb)
utr3 <- threeUTRsByTranscript(txdb)
introns <- intronsByTranscript(txdb)

# Annotate the BED regions
exon_hits <- findOverlaps(bed_gr, exons)
bed_gr$annotation <- "intergenic"
bed_gr$annotation[queryHits(exon_hits)] <- "exon"

utr5_hits <- findOverlaps(bed_gr, utr5)
bed_gr$annotation[queryHits(utr5_hits)] <- "5'UTR"

utr3_hits <- findOverlaps(bed_gr, utr3)
bed_gr$annotation[queryHits(utr3_hits)] <- "3'UTR"

intron_hits <- findOverlaps(bed_gr, introns)
bed_gr$annotation[queryHits(intron_hits)] <- "intron"

# Convert to data frame for output
bed_df <- as.data.frame(bed_gr)
# Convert BED start positions to one-based by adding 1
bed_df$start <- bed_df$start - 1

# Ensure 'annotation' is included in the output
output_bed <- file.path(args$output, paste0(args$name, "_annotated.bed"))
write.table(bed_df[, c("seqnames", "start", "end", "annotation")], file = output_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Load the TSV file
tsv_df <- read.table(args$tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Merge the BED and TSV files based on the 'chrom', 'start', and 'end' columns
merged_df <- merge(bed_df, tsv_df, by.x = c("seqnames", "start", "end"), by.y = c("chrom", "start", "end"), all.x = TRUE)

# Save the merged output as TSV
output_merged <- file.path(args$output, paste0(args$name, "_annotated.tsv"))
write.table(merged_df, file = output_merged, sep = "\t", quote = FALSE, row.names = FALSE)

# Prepare data for plotting
annotation_counts <- table(bed_df$annotation)
annotation_df <- as.data.frame(annotation_counts)
colnames(annotation_df) <- c("Feature", "Count")

# Create the annotation plot
plot <- ggplot(annotation_df, aes(x = Feature, y = Count, fill = Feature)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  labs(title = "Distribution of Genomic Features",
       x = "Genomic Feature",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = Count), vjust = -0.5, size = 3.5)

# Save the plot
output_plot1 <- file.path(args$output, paste0(args$name, "_genomic_features_distribution_barplot.png"))
ggsave(output_plot1, plot = plot, width = 8, height = 6)

# Load input for second script
queryFile <- args$bed
outputFile <- file.path(args$output, paste0(args$name, "_genomic_and_chr_distribution.pdf"))

query <- import(queryFile)

# Plot chromosome distribution
x <- calcChromBinsRef(query, "hg38")
chrplot <- plotChromBins(x)

# Calculate genomic partition distribution and plot
gp <- calcPartitionsRef(query, "hg38")
q <- ggplot(gp, aes(x = partition, y = Freq, fill = partition)) +
  geom_bar(stat = "identity") +
  ggtitle("Genomic features distribution plot") +
  rremove("legend") + 
  rremove("xlab")

# Expected partition distribution plots
ep <- calcExpectedPartitionsRef(query, "hg38")
p <- plotExpectedPartitions(ep) +
  aes(fill = partition) +
  ggtitle("Expected genomic distribution plots") +
  labs(subtitle = "Distribution of intergenic vs promoter space 
in the genome is not uniform. 
Here, we produce plots showing the log10(Obs/Exp) of the 
regions across genomic partitions.") +
  rremove("legend") +
  theme(
    plot.title = element_text(color = "red", size = 14, face = "bold"),
    plot.subtitle = element_text(color = "blue"),
    plot.caption = element_text(color = "green", face = "italic")
  )

# Combine plots
myplot <- ggarrange(chrplot, 
                    ggarrange(q, p, ncol = 2, labels = c("B", "C")),
                    nrow = 2, 
                    labels = "A")

# Export combined plot
ggpubr::ggexport(myplot, filename = outputFile, width = 14, height = 12)

cat("Annotated BED file saved to:", output_bed, "\n")
cat("Annotated TSV file saved to:", output_merged, "\n")
cat("Genomic features distribution plot saved to:", output_plot1, "\n")
cat("Genomic distribution plots saved to:", outputFile, "\n")

