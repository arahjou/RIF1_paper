# Install required packages if not already installed
if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("GenomicRanges")
}

if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}

# Load necessary libraries
library(GenomicRanges)
library(VennDiagram)



# Function to read a BED file and convert it to a GRanges object
read_bed <- function(file) {
  df <- read.table(file, header = FALSE)
  gr <- GRanges(seqnames = df[, 1],
                ranges = IRanges(start = df[, 2], end = df[, 3]))
  return(gr)
}


# Read the two BED files
bed1 <- read_bed("RIF1_paper_20240901/AID off-target/AID_Offtarget_mm10.bed")
bed2 <- read_bed("RIF1_paper_20240901/AID off-target/RIF1_new_peaks.bed")
bed3 <- read_bed("RIF1_paper_20240901/AID off-target/BLIP1_peaks.bed")


# Find overlaps
overlaps <- findOverlaps(bed1, bed2)

# Identify the number of unique peaks in each file and the overlaps
unique_peaks_file1 <- length(bed1)
unique_peaks_file2 <- length(bed2)
num_overlaps <- length(unique(queryHits(overlaps)))

# Non-overlapping counts
num_non_overlaps_1 <- unique_peaks_file1 - num_overlaps
num_non_overlaps_2 <- unique_peaks_file2 - length(unique(subjectHits(overlaps)))

# Create a Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = unique_peaks_file1,
  area2 = unique_peaks_file2,
  cross.area = num_overlaps,
  category = c("AID-Offtarget", "RIF1_peaks"),
  fill = c("dodgerblue", "tomato"),
  alpha = 0.5,
  cat.pos = c(-20, 20),
  cat.dist = 0.03,
  lty = "blank",
  cex = 1,
  cat.cex = 1,
  cat.col = c("dodgerblue", "tomato")
)

# Save the Venn diagram as a file
pdf("/Users/probio/Desktop/RIF1_paper_20240901/AID off-target/venn_diagram.pdf", width=7, height=7)
grid.draw(venn.plot)
dev.off()

# 
# --------------------------------------------


# Create a GRanges object representing the overlapping regions between bed1 and bed2
overlap_granges <- pintersect(bed1[queryHits(overlaps)], bed2[subjectHits(overlaps)])
overlap_granges <- reduce(overlap_granges) # Merge overlapping ranges into a single range

# Find overlaps between bed3 and the overlapping regions
overlaps_with_bed3 <- findOverlaps(bed3, overlap_granges)

# Identify the number of unique peaks in bed3 and the overlaps with the overlap_granges
unique_peaks_bed3 <- length(bed3)
num_overlaps_with_bed3 <- length(unique(queryHits(overlaps_with_bed3)))

# Non-overlapping counts for bed3 and overlap_granges
num_non_overlaps_bed3 <- unique_peaks_bed3 - num_overlaps_with_bed3
num_non_overlaps_overlap_granges <- length(overlap_granges) - num_overlaps_with_bed3

# Create a Venn diagram for bed3 and overlap_granges
venn.plot2 <- draw.pairwise.venn(
  area1 = length(overlap_granges),
  area2 = unique_peaks_bed3,
  cross.area = num_overlaps_with_bed3,
  category = c("Overlap from AID-off and RIF1", "Blimp1_peaks"),
  fill = c("purple", "orange"),
  alpha = 0.5,
  cat.pos = c(-20, 20),
  cat.dist = 0.03,
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("purple", "orange")
)

# Save the Venn diagram as a file
pdf("/Users/probio/Desktop/RIF1_paper_20240901/AID off-target/venn_diagram_Blimp1.pdf")
grid.draw(venn.plot2)
dev.off()

#--------

# Install required packages if not already installed
if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("GenomicRanges")
}

if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}

# Load necessary libraries
library(GenomicRanges)
library(VennDiagram)

# Function to read a BED file and convert it to a GRanges object
read_bed <- function(file) {
  df <- read.table(file, header = FALSE)
  gr <- GRanges(seqnames = df[, 1],
                ranges = IRanges(start = df[, 2], end = df[, 3]))
  return(gr)
}

# Read the three BED files
bed1 <- read_bed("/Users/arahjou/Downloads/RIF1_paper_20240901/AID off-target/AID_Offtarget_mm10.bed")
bed2 <- read_bed("/Users/arahjou/Downloads/RIF1_paper_20240901/AID off-target/RIF1_new_peaks.bed")
bed3 <- read_bed("/Users/arahjou/Downloads/RIF1_paper_20240901/AID off-target/BLIP1_peaks.bed")

# Find overlaps between the three BED files
overlap1_2 <- intersect(bed1, bed2)
overlap1_3 <- intersect(bed1, bed3)
overlap2_3 <- intersect(bed2, bed3)
overlap_all <- Reduce(intersect, list(bed1, bed2, bed3))

# Save the overlapping regions to a file
write.table(as.data.frame(overlap_all), "/Users/arahjou/Downloads/RIF1_paper_20240901/AID off-target/overlap_regions.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# Plot a Venn diagram
venn <- venn.diagram(
  x = list(
    AID = bed1,
    RIF1 = bed2,
    BLIP1 = bed3
  ),
  category.names = c("AID", "RIF1", "BLIP1"),
  filename = NULL, # NULL to plot directly
  output = TRUE
)

# Save Venn diagram
png("overlap_venn_diagram.png")
grid.draw(venn)
dev.off()

# Display the number of overlaps
cat("Number of overlaps between AID and RIF1:", length(overlap1_2), "\n")
cat("Number of overlaps between AID and BLIP1:", length(overlap1_3), "\n")
cat("Number of overlaps between RIF1 and BLIP1:", length(overlap2_3), "\n")
cat("Number of overlaps between all three:", length(overlap_all), "\n")
