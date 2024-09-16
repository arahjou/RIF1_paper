# Load necessary libraries
library(biomaRt)
library(Gviz)
library(rtracklayer)

# Set up the Ensembl biomaRt object
mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

# Get the gene coordinates for the Aicda gene
gene_info <- getBM(
  attributes = c('chromosome_name', 'start_position', 'end_position', 'external_gene_name'),
  filters = 'external_gene_name', values = 'Aicda', mart = mart
)

# Extract the start, end, and chromosome information for the Aicda gene
start_pos <- gene_info$start_position[1]
end_pos <- gene_info$end_position[1]
chromosome <- gene_info$chromosome_name[1]  # Likely "6"

# Create the GeneRegionTrack for Aicda gene
gene_track <- BiomartGeneRegionTrack(
  genome = "mm10",            # UCSC genome identifier
  chromosome = chromosome,    # Chromosome name from gene_info
  start = start_pos,          # Start position from gene_info
  end = end_pos,              # End position from gene_info
  biomart = mart,             # Use the manually defined biomaRt object
  name = "Aicda Gene",
  showId = TRUE,
  transcriptAnnotation = "transcript"
)

# Load the .bw files
bw1 <- import("/Users/arahjou/Downloads/RIF1_paper_20240901/Peak_visualization/RIF1_ChIP.bigwig", format = "BigWig")
bw2 <- import("/Users/arahjou/Downloads/RIF1_paper_20240901/Peak_visualization/Mock.bigwig", format = "BigWig")
#bw3 <- import("file3.bw", format = "BigWig")
#bw4 <- import("file4.bw", format = "BigWig")

# Load the .bed files
bed1 <- import("/Users/arahjou/Downloads/RIF1_paper_20240901/Peak_visualization/RIF1_Peaks.bed", format = "BED")
#bed2 <- import("file2.bed", format = "BED")

# Create DataTracks for the bigWig files
track_bw1 <- DataTrack(bw1, name = "Track 1", type = "l", genome = "mm10")
track_bw2 <- DataTrack(bw2, name = "Track 2", type = "l", genome = "mm10")
#track_bw3 <- DataTrack(bw3, name = "Track 3", type = "l", genome = "mm10")
#track_bw4 <- DataTrack(bw4, name = "Track 4", type = "l", genome = "mm10")

# Create AnnotationTracks for the bed files
track_bed1 <- AnnotationTrack(bed1, name = "BED 1", genome = "mm10")
#track_bed2 <- AnnotationTrack(bed2, name = "BED 2", genome = "mm10")

# Plot the gene and tracks
plotTracks(list(
  GenomeAxisTrack(),
  gene_track,            # Include the gene track for Aicda
  track_bw1, track_bw2, #track_bw3, track_bw4,
  track_bed1, #track_bed2
), from = start_pos, to = end_pos, chromosome = chromosome)


# Plot the gene and tracks without empty arguments
plotTracks(list(
  GenomeAxisTrack(),
  gene_track,            # Include the gene track for Aicda
  track_bw1, track_bw2,  # Only two bigWig tracks included
  track_bed1             # Only one BED track included
), from = start_pos, to = end_pos, chromosome = chromosome)
??plotTracks


#----------------------------

# Load necessary libraries
library(biomaRt)
library(Gviz)
library(rtracklayer)

# Set up the Ensembl biomaRt object
mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

# Get the gene coordinates for the Aicda gene
gene_info <- getBM(
  attributes = c('chromosome_name', 'start_position', 'end_position', 'external_gene_name'),
  filters = 'external_gene_name', values = 'Aicda', mart = mart
)

# Extract the start, end, and chromosome information for the Aicda gene
start_pos <- gene_info$start_position[1]
end_pos <- gene_info$end_position[1]
chromosome <- gene_info$chromosome_name[1]  # Likely "6"

# Define a specific region around the Aicda gene, e.g., +/- 5000 bases around the gene
region_start <- max(1, start_pos - 5000)  # Ensuring start position is not negative
region_end <- end_pos + 5000              # Add 5000 bases downstream of the gene

# Create the GeneRegionTrack for Aicda gene
gene_track <- BiomartGeneRegionTrack(
  genome = "mm10",            # UCSC genome identifier
  chromosome = chromosome,    # Chromosome name from gene_info
  start = start_pos,          # Start position from gene_info
  end = end_pos,              # End position from gene_info
  biomart = mart,             # Use the manually defined biomaRt object
  name = "Aicda Gene",
  showId = TRUE,
  transcriptAnnotation = "transcript"
)

# Load the .bw files
bw1 <- import("/Users/arahjou/Downloads/RIF1_paper_20240901/Peak_visualization/RIF1_ChIP.bigwig", format = "BigWig")
bw2 <- import("/Users/arahjou/Downloads/RIF1_paper_20240901/Peak_visualization/Mock.bigwig", format = "BigWig")

# Load the .bed files
bed1 <- import("/Users/arahjou/Downloads/RIF1_paper_20240901/Peak_visualization/RIF1_Peaks.bed", format = "BED")

# Create DataTracks for the bigWig files
track_bw1 <- DataTrack(bw1, name = "Track 1", type = "l", genome = "mm10")
track_bw2 <- DataTrack(bw2, name = "Track 2", type = "l", genome = "mm10")

# Create AnnotationTracks for the bed files
track_bed1 <- AnnotationTrack(bed1, name = "BED 1", genome = "mm10")

# Plot the gene and tracks for the specific region
plotTracks(list(
  GenomeAxisTrack(),
  gene_track,            # Include the gene track for Aicda
  track_bw1, track_bw2,  # Only two bigWig tracks included
  track_bed1             # Only one BED track included
), from = region_start, to = region_end, chromosome = chromosome)


# -------------------

# Load necessary libraries
library(biomaRt)
library(Gviz)
library(rtracklayer)

# Set up the Ensembl biomaRt object
mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

# Get the gene coordinates for the Aicda gene
gene_info <- getBM(
  attributes = c('chromosome_name', 'start_position', 'end_position', 'external_gene_name'),
  filters = 'external_gene_name', values = 'Aicda', mart = mart
)

# Extract the start, end, and chromosome information for the Aicda gene
start_pos <- gene_info$start_position[1]
end_pos <- gene_info$end_position[1]
chromosome <- gene_info$chromosome_name[1]  # Likely "6"

# Define a specific region around the Aicda gene, e.g., +/- 5000 bases around the gene
region_start <- max(1, start_pos - 5000)  # Ensuring start position is not negative
region_end <- end_pos + 5000              # Add 5000 bases downstream of the gene

# Create the GeneRegionTrack for Aicda gene
gene_track <- BiomartGeneRegionTrack(
  genome = "mm10",            # UCSC genome identifier
  chromosome = chromosome,    # Chromosome name from gene_info
  start = start_pos,          # Start position from gene_info
  end = end_pos,              # End position from gene_info
  biomart = mart,             # Use the manually defined biomaRt object
  name = "Aicda Gene",
  showId = TRUE,
  transcriptAnnotation = "transcript"
)

# Load the .bw files
bw1 <- import("/Users/arahjou/Downloads/RIF1_paper_20240901/Peak_visualization/RIF1_ChIP.bigwig", format = "BigWig")
bw2 <- import("/Users/arahjou/Downloads/RIF1_paper_20240901/Peak_visualization/Mock.bigwig", format = "BigWig")

# Load the .bed files
bed1 <- import("/Users/arahjou/Downloads/RIF1_paper_20240901/Peak_visualization/RIF1_Peaks.bed", format = "BED")

# Create DataTracks for the bigWig files, with fixed ylim (0 to 50)
track_bw1 <- DataTrack(bw1, name = "Track 1", type = "polygon", genome = "mm10", ylim = c(0, 50))
track_bw2 <- DataTrack(bw2, name = "Track 2", type = "polygon", genome = "mm10", ylim = c(0, 50))

# Create AnnotationTracks for the bed files
track_bed1 <- AnnotationTrack(bed1, name = "BED 1", genome = "mm10")

# Plot the gene and tracks for the specific region
plotTracks(list(
  GenomeAxisTrack(),
  gene_track,            # Include the gene track for Aicda
  track_bw1, track_bw2,  # Only two bigWig tracks included, both with fixed y-axis height
  track_bed1             # Only one BED track included
), from = region_start, to = region_end, chromosome = chromosome)

??DataTrack
# -------------------

# Load necessary libraries
library(biomaRt)
library(Gviz)
library(rtracklayer)

# Set up the Ensembl biomaRt object
mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

# Get the gene coordinates for the Aicda gene
gene_info <- getBM(
  attributes = c('chromosome_name', 'start_position', 'end_position', 'external_gene_name'),
  filters = 'external_gene_name', values = 'Exosc3', mart = mart
)

# Extract the start, end, and chromosome information for the Aicda gene
start_pos <- gene_info$start_position[1]
end_pos <- gene_info$end_position[1]
chromosome <- gene_info$chromosome_name[1]  # Likely "6"

# Define a specific region around the Aicda gene, e.g., +/- 5000 bases around the gene
region_start <- max(1, start_pos - 10000)  # Ensuring start position is not negative
region_end <- end_pos + 20000              # Add 5000 bases downstream of the gene

# Create the GeneRegionTrack for Aicda gene
gene_track <- BiomartGeneRegionTrack(
  genome = "mm10",            # UCSC genome identifier
  chromosome = chromosome,    # Chromosome name from gene_info
  start = start_pos,          # Start position from gene_info
  end = end_pos,              # End position from gene_info
  biomart = mart,             # Use the manually defined biomaRt object
  name = "Exosc3 gene Gene",
  showId = TRUE,
  transcriptAnnotation = "transcript"
)

# Load the .bw files
bw1 <- import("/Users/arahjou/Downloads/RIF1_paper_20240901/Peak_visualization/RIF1_ChIP.bigwig", format = "BigWig")
bw2 <- import("/Users/arahjou/Downloads/RIF1_paper_20240901/Peak_visualization/Mock.bigwig", format = "BigWig")

# Load the .bed files
bed1 <- import("/Users/arahjou/Downloads/RIF1_paper_20240901/Peak_visualization/RIF1_Peaks.bed", format = "BED")

# Create DataTracks for the bigWig files, with filled areas under the lines
# Custom line and fill colors using `col` for border and `col.fill` for the fill
track_bw1 <- DataTrack(bw1, name = "RIF1", type = "hist", genome = "mm10", ylim = c(0, 50), 
                       col.histogram = "darkblue", fill = "darkblue")  # Dark blue line, light blue fill
track_bw2 <- DataTrack(bw2, name = "Mock", type = "hist", genome = "mm10", ylim = c(0, 50), 
                       col.histogram ="red" , fill = "red")        # Dark red line, pink fill

# Create AnnotationTracks for the bed files, with custom colors
track_bed1 <- AnnotationTrack(bed1, name = "Peaks", genome = "mm10", 
                              fill = "lightgreen", col = "darkgreen")  # Light green fill, dark green border

# Plot the gene and tracks for the specific region
plotTracks(list(
  GenomeAxisTrack(),
  gene_track,            # Include the gene track for Aicda
  track_bw1, track_bw2,  # Filled area under the BigWig tracks with custom colors
  track_bed1             # BED track with custom colors
), from = region_start, to = region_end, chromosome = chromosome)



# ------------ This version works with just gene location on the genome instead of gene name.

# Define a specific region: chr4:45316613-45320629
chromosome <- "12"          # Define the chromosome of interest (chr4)
region_start <- 113314978    # Define the start position
region_end <- 113328929   # Define the end position

#chr12:113421971-113428142 mm10_Smu
#chr12:113314978-113328929 mm10_Sg
#chr10:44394100-44567131
#chr4:45316455-45329404
# Create the GeneRegionTrack for the specified region
gene_track <- BiomartGeneRegionTrack(
  genome = "mm10",            # UCSC genome identifier
  chromosome = chromosome,    # Chromosome name ("4" in this case)
  start = region_start,       # Start position for the region
  end = region_end,           # End position for the region
  biomart = mart,             # Use the manually defined biomaRt object
  name = "Ighm",       # Change the gene name if needed
  showId = TRUE,
  transcriptAnnotation = "transcript"
)

# Load the .bw files
bw1 <- import("/Users/arahjou/Downloads/RIF1_paper_20240901/Peak_visualization/RIF1_ChIP.bigwig", format = "BigWig")
bw2 <- import("/Users/arahjou/Downloads/RIF1_paper_20240901/Peak_visualization/Mock.bigwig", format = "BigWig")

# Load the .bed files
#bed1 <- import("/Users/arahjou/Downloads/RIF1_paper_20240901/Peak_visualization/RIF1_Peaks.bed", format = "BED")
bed1 <- import("/Users/arahjou/Downloads/RIF1_paper_20240901/Peak_visualization/20190820_IgHmm10.bed", format = "BED")

# Create DataTracks for the bigWig files, with filled areas under the lines
track_bw1 <- DataTrack(bw1, name = "RIF1", type = "hist", genome = "mm10", ylim = c(0, 60), 
                       col.histogram = "darkblue", fill = "darkblue")
track_bw2 <- DataTrack(bw2, name = "Mock", type = "hist", genome = "mm10", ylim = c(0, 60), 
                       col.histogram ="red", fill = "red")

# Create AnnotationTracks for the bed files, with custom colors
track_bed1 <- AnnotationTrack(bed1, name = "Peaks", genome = "mm10", 
                              fill = "lightgreen", col = "darkgreen")

# Plot the gene and tracks for the specified region
plotTracks(list(
  GenomeAxisTrack(),
  gene_track,            # Include the gene track for the specified region
  track_bw1, track_bw2,  # Filled area under the BigWig tracks with custom colors
  track_bed1             # BED track with custom colors
), from = region_start, to = region_end, chromosome = chromosome)




#---------------- Plotting size of RIF1 peaks


Blimp1 <-  read.table("/Users/arahjou/Downloads/RIF1_paper_20240901/RIF1_two/Peaks/BLIP1_peaks.txt", header = FALSE)
Blimp1_peak_size <- Blimp1$V3 - Blimp1$V2
summary(Blimp1_peak_size)
Rif1 <- read.table("/Users/arahjou/Downloads/RIF1_paper_20240901/RIF1_two/Peaks/RIF1_peaks.txt", header = FALSE)
Rif1_peak_size <- Rif1$V3 - Rif1$V2
summary(Rif1_peak_size)

hist(Blimp1_peak_size, breaks = 100, col = "lightblue", main = "Distribution of Peak Sizes (Blimp1)", xlab = "Peak Size")
hist(Rif1_peak_size, breaks = 100, col = "lightgreen", main = "Distribution of Peak Sizes (Rif1)", xlab = "Peak Size")

# join the two histograms
hist(Blimp1_peak_size, breaks = 100, col = "lightblue", main = "Distribution of Peak Sizes", xlab = "Peak Size")
hist(Rif1_peak_size, breaks = 100, col = "lightgreen", add = TRUE)


# create a long dataframe from the two peak size vectors

peak_sizes <- data.frame(
  PeakSize = c(Blimp1_peak_size, Rif1_peak_size),
  Factor = factor(rep(c("Blimp1", "Rif1"), c(length(Blimp1_peak_size), length(Rif1_peak_size))))
)

# Create a boxplot of the peak sizes using ggplot2
library(ggplot2)
ggplot(peak_sizes, aes(x = Factor, y = PeakSize, fill = Factor)) +
  geom_boxplot() +
  labs(title = "Distribution of Peak Sizes", x = "Factor", y = "Peak Size") +
  theme_minimal() +
  scale_fill_manual(values = c("lightblue", "lightgreen"))

# violin plot ggplo2
ggplot(peak_sizes, aes(x = Factor, y = PeakSize, fill = Factor)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "Distribution of Peak Sizes", x = "Factor", y = "Peak Size") +
  ylim(0, 4000) +
  #y axis label to peak size (bp)
  ylab("Peak Size (bp)") +
  theme_minimal() +
  scale_fill_manual(values = c("lightblue", "lightgreen"))
ggsave("/Users/arahjou/Downloads/RIF1_paper_20240901/RIF1_two/Peaks/peak_sizes_violin_plot.pdf", width = 6, height = 7, dpi = 300)

library(dplyr)
peak_sizes$PeakSize
peak_sizes %>% 
  group_by(Factor) %>% 
  summarize(mean_size = mean(PeakSize, na.rm = TRUE),
            median_size = median(PeakSize, na.rm = TRUE),
            sd_size = sd(PeakSize, na.rm = TRUE),
            min_size = min(PeakSize, na.rm = TRUE),
            max_size = max(PeakSize, na.rm = TRUE))



  # making IgH image--------------------
# Load necessary libraries
library(biomaRt)
library(Gviz)
library(rtracklayer)

# Set up the Ensembl biomaRt object
mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

# Define a specific region: chr4:45316613-45320629
chromosome <- "12"          # Define the chromosome of interest (chr4)
region_start <- 114552337    # Define the start position
region_end <- 114565934   # Define the end position

#chr12:114658069-114667335 # mm9 Smu
#chr12:114552337-114565934 # mm9 Sg
#chr10:44394100-44567131
#chr4:45316455-45329404
# Create the GeneRegionTrack for the specified region
gene_track <- BiomartGeneRegionTrack(
  genome = "mm9",            # UCSC genome identifier
  chromosome = chromosome,    # Chromosome name ("4" in this case)
  start = region_start,       # Start position for the region
  end = region_end,           # End position for the region
  biomart = mart,             # Use the manually defined biomaRt object
  name = "Ighm",       # Change the gene name if needed
  showId = TRUE,
  transcriptAnnotation = "transcript"
)

# Load the .bw files
bw1 <- import("/Users/arahjou/Desktop/Robert_rif1_mm9/rif1_pooled.bw", format = "BigWig")
bw2 <- import("/Users/arahjou/Desktop/Robert_rif1_mm9/rif1_mockup.bw", format = "BigWig")

# Load the .bed files
bed1 <- import("/Users/arahjou/Desktop/Robert_rif1_mm9/IgH_features_mm9.bed", format = "BED")

# Create DataTracks for the bigWig files, with filled areas under the lines
track_bw1 <- DataTrack(bw1, name = "RIF1", type = "hist", genome = "mm9", ylim = c(0, 2), 
                       col.histogram = "darkblue", fill = "darkblue")
track_bw2 <- DataTrack(bw2, name = "Mock", type = "hist", genome = "mm9", ylim = c(0, 2), 
                       col.histogram ="red", fill = "red")

# Create AnnotationTracks for the bed files, with custom colors
track_bed1 <- AnnotationTrack(bed1, name = "Peaks", genome = "mm9", 
                              fill = "lightgreen", col = "darkgreen")

# Plot the gene and tracks for the specified region
plotTracks(list(
  GenomeAxisTrack(),
  gene_track,            # Include the gene track for the specified region
  track_bw1, track_bw2,  # Filled area under the BigWig tracks with custom colors
  track_bed1             # BED track with custom colors
), from = region_start, to = region_end, chromosome = chromosome)




