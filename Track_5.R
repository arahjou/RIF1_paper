####################Gene track#################


library(rtracklayer)
library(GenomicRanges)
library(dplyr)

# Load the GTF file - replace with the path to your GTF file
gtf_path <- "/Users/probio/Downloads/gencode.vM10.annotation.gtf"
gtf_data <- import(gtf_path)

# Filter for Nedd4 gene annotations
# Replace "Nedd4" with the exact gene name or ID used in your GTF file
nedd4_gene <- subset(gtf_data, gene_name == "Sell")


# Extract exons for Nedd4
nedd4_exons <- subset(nedd4_gene, type == "exon")

# Create a data frame for exons
nedd4_exons_df <- data.frame(
  start = start(nedd4_exons),
  end = end(nedd4_exons)
)

# Get intron regions by finding gaps between exons
nedd4_introns <- GenomicRanges::gaps(nedd4_exons)

# Create a data frame for introns
nedd4_introns_df <- data.frame(
  start = start(nedd4_introns),
  end = end(nedd4_introns)
)


# Filter for exons that belong to the specific transcript ID
specific_transcript_exons <- subset(nedd4_gene, mcols(nedd4_gene)$transcript_id == "ENSMUST00000192047.5")

# Create a data frame for exons of the specific transcript
specific_transcript_exons_df <- data.frame(
  start = start(specific_transcript_exons),
  end = end(specific_transcript_exons),
  strand = strand(specific_transcript_exons),
  score = mcols(specific_transcript_exons)$score
)

# You can also check the number of exons for this specific transcript
num_exons_specific_transcript <- length(specific_transcript_exons)

# Print the number of exons to the console
print(num_exons_specific_transcript)



####################bigwig_file#################

# Load necessary libraries
library(Signac)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(cowplot) # Optional, for arranging multiple plots

# Define file paths
blimp1_bw_path <- "/Users/probio/Downloads/BLIMP1.bigwig"
rif1_bw_path <- "/Users/probio/Downloads/RIF1.bigwig"

# Load BigWig files
blimp1_bw <- import(blimp1_bw_path)
rif1_bw <- import(rif1_bw_path)



# Assuming you know the genomic region of interest for 'Nedd4' gene in mm10
# Example: chrX:10,000,000-11,000,000
# Please replace with the actual coordinates of the Nedd4 gene
nedd4_region <- GRanges("chr1", IRanges(start=164059982, end=164086181))

# Subset the tracks to the Nedd4 region
blimp1_bw_nedd4 <- subsetByOverlaps(blimp1_bw, nedd4_region)
rif1_bw_nedd4 <- subsetByOverlaps(rif1_bw, nedd4_region)


# Assuming you've already prepared the data frames for each of the tracks (blimp1_bw_df, rif1_bw_df, blimp1_peaks_df, rif1_peaks_df)
# and you have a separate data frame for Nedd4 exons (nedd4_exons_df) and introns (nedd4_introns_df)

blimp1_bw_df = as.data.frame(blimp1_bw_nedd4)
rif1_bw_df = as.data.frame(rif1_bw_nedd4)


# Plot for BLIMP1 signal with filled area under the peaks
p1 <- ggplot(blimp1_bw_df, aes(x = start, y = score)) +
  geom_area(fill = "#9a11d9", alpha = 0.5) +  # Fills under the line
  geom_line(color = "#9a11d9") +  # Keeps the line on top of the fill
  labs(title = "BLIMP1 Signal") +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())

# Plot for RIF1 signal with filled area under the peaks
p2 <- ggplot(rif1_bw_df, aes(x = start, y = score)) +
  geom_area(fill = "#ff792b", alpha = 0.5) +  # Fills under the line
  geom_line(color = "#ff792b") +  # Keeps the line on top of the fill
  labs(title = "RIF1 Signal") +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())


# Create a data frame for exons of the specific transcript
specific_transcript_exons_df <- data.frame(
  start = start(specific_transcript_exons),
  end = end(specific_transcript_exons),
  strand = strand(specific_transcript_exons),
  score = mcols(specific_transcript_exons)$score
)


specific_transcript_exons_df = distinct(specific_transcript_exons_df)
specific_transcript_exons_df <- specific_transcript_exons_df[seq_along(specific_transcript_exons_df[ ,1]) %% 2 == 0, ]

#summary(specific_transcript_exons_df)
# Plot for Nedd4 gene (exons and introns)
# Assuming nedd4_exons_df and nedd4_introns_df are prepared



p3 <- ggplot() +
  geom_segment(data=specific_transcript_exons_df, aes(x=start, xend=end, y=0, yend=0), color="red", linewidth=3) +
  labs(title = "Nedd4 Gene (Exons and Introns)") +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank()) +
  ylim(-1, 1) +
  geom_hline(yintercept=0, size=0.5, color = 'gray') # Adding horizontal line at y=0 with size 1

# You can now plot p3 or save it to a file


# Arrange all plots
plot_grid(p1, p2, p3, ncol=1, align='v')



##################
# Plot for BLIMP1 signal with filled area under the peaks
p1 <- ggplot(blimp1_bw_df, aes(x = start, y = score)) +
  geom_area(fill = "black", alpha = 0.5) +  # Fills under the line
  geom_line(color = "black") +  # Keeps the line on top of the fill
  labs(title = "BLIMP1 Signal") +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())

# Plot for RIF1 signal with filled area under the peaks
p2 <- ggplot(rif1_bw_df, aes(x = start, y = score)) +
  geom_area(fill = "black", alpha = 0.5) +  # Fills under the line
  geom_line(color = "black") +  # Keeps the line on top of the fill
  labs(title = "RIF1 Signal") +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank())




p3 <- ggplot() +
  geom_segment(data=Nedd4_exons, aes(x=start, xend=end, y=0, yend=0), color="black", linewidth=3) +
  labs(title = "Nedd4 Gene (Exons and Introns)") +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank()) +
  ylim(-1, 1) +
  geom_hline(yintercept=0, size=0.5, color = 'black') # Adding horizontal line at y=0 with size 1

# You can now plot p3 or save it to a file


# Arrange all plots
plot_grid(p1, p2, p3, ncol=1, align='v')

