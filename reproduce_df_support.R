# Reproduce Dataframe Support for plot_taxa_star

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(phyloseq)
library(microViz)

# Source the modified function
source("R/utils-radar.R")
source("R/utils-colors.R")
source("R/taxa_star_plot.R")

# 1. Create Mock Data
set.seed(123)
n_samples <- 10
taxa_names <- c("Firmicutes", "Bacteroidetes", "Proteobacteria", "Actinobacteria")
other_taxa <- c("T5", "T6")

# Create abundance data (counts)
abundance_matrix <- matrix(sample(1:100, n_samples * (length(taxa_names) + length(other_taxa)), replace = TRUE),
                           nrow = n_samples,
                           dimnames = list(paste0("Sample", 1:n_samples), c(taxa_names, other_taxa)))

# Create metadata
metadata <- data.frame(
  SampleID = paste0("Sample", 1:n_samples),
  Group = sample(c("Control", "Treatment"), n_samples, replace = TRUE),
  Location = sample(c("Gut", "Skin"), n_samples, replace = TRUE)
)

# Combine into a single dataframe
df <- cbind(metadata, as.data.frame(abundance_matrix))

print("Dataframe created:")
print(head(df))

# 2. Test plot_taxa_star with Dataframe
tryCatch({
  p <- plot_taxa_star(
    physeq = df,
    sample_var = "Group",
    taxa_rank = "Phylum", # This will be the name of the rank column in the dummy tax table
    taxa_names = taxa_names[1:3],
    samplecolumn = "SampleID"
  )
  
  print("Plot generated successfully!")
  
  # Test with SE and custom alpha
  p_se <- plot_taxa_star(
    physeq = df,
    sample_var = "Group",
    taxa_rank = "Phylum",
    taxa_names = taxa_names[1:3],
    samplecolumn = "SampleID",
    error_bar = "SE",
    fill_alpha = 0.2
  )
  print("Plot (SE, Alpha 0.2) generated successfully!")
  
  # Test Default Taxa Selection (Top 4)
  p_default <- plot_taxa_star(
    physeq = df,
    sample_var = "Group",
    taxa_rank = "Phylum",
    # taxa_names argument omitted (should default to top 4)
    samplecolumn = "SampleID"
  )
  print("Plot (Default Top 4 Taxa) generated successfully!")
  
  print(class(p))
  
  # Check if it is a ggplot object
  if (inherits(p, "ggplot")) {
    print("Verification PASSED: Result is a ggplot object.")
  } else {
    print("Verification FAILED: Result is not a ggplot object.")
  }
  
}, error = function(e) {
  print(paste("Verification FAILED with error:", e$message))
  traceback()
})
