# microbiomestar

<!-- badges: start -->
<!-- badges: end -->

The `microbiomestar` package provides unique visualization techniques for microbiome data using `phyloseq` and `microViz`. Its flagship features involve replacing standard points or bars with "star" (radar) plots that show the relative abundance of core taxa or any taxa of interest. These multi-dimensional shapes can be projected onto PCoA plots or arranged in core matrices to easily visualize differences in community structures across samples or treatments (up to 4 groups/treatments for comparative core matrices).

<img width="1315" height="747" alt="image" src="https://github.com/user-attachments/assets/f33bea51-a39e-47d6-9680-9a998c4c031b" />
*Enterotypes plotted using data from cvanmf (plot_taxa_star())*
<img width="1315" height="747" alt="image" src="https://github.com/user-attachments/assets/89dd71a7-bb12-43b5-979c-83a2a9a43128" />
*PCoA axes plotted using data from the GlobalPatterns dataset (plot_pcoa_star())*


## Installation

You can install the development version of `microbiomestar` from GitHub with:
``` r
BiocManager::install("kadenmuffett/microbiomestarplots")
```

## Features

### `plot_taxa_star()`
Creates radar/star plots to visualize the relative abundance of specified taxa (or the top N highest-abundant taxa automatically) across different sample groups, complete with standard error or IQR shading overlays.

### `plot_pcoa_star()`
Performs Principal Coordinate Analysis (PCoA) on beta-diversity distance matrices (such as Bray-Curtis) and replaces standard geometric points on the scatterplot with star plots representing the taxonomic breakdown of each sample at its ordinal coordinate.

### `plot_core_matrix()`
Generates a facet matrix displaying the "unique core" ASVs across groups. Shows the abundance oftaxa that are considered "core" to one treatment/group but completely missing from the core of the other groups.

### `coreselect()`
Provides a presence/absence comparison tool for assessing the core microbiome thresholds across various groupings within your phyloseq object. This tool is exploratory and can be used for selection of appropriate taxa for plot_taxa_star() or to assess the appropriateness of plot_core_matrix() for your data.
### *Please access the vignette to test it out, as it includes data appropriate for plot_core_matrix and coreselect, GlobalPatterns is not.*
## Examples

```r
library(microbiomestar)
library(phyloseq)

data("GlobalPatterns")
# Subset to a smaller physeq for demonstration
gp <- subset_samples(GlobalPatterns, SampleType %in% c("Feces", "Ocean", "Soil"))
# For better and more complete examples, access the vignettes
# vignette("microbiomeplots")
# Plot Taxa Star Plot with Standard Error shading
plot_taxa_star(
  physeq = gp,
  sample_var = "SampleType",
  taxa_rank = "Phylum",
  samplecolumn = "X.SampleID",
  error_bar = "SE"
)

# Plot PCoA mapped with Star Plots based on Bray-Curtis distance!
plot_pcoa_star(
  physeq = gp,
  sample_var = "SampleType",
  taxa_rank = "Phylum",
  samplecolumn = "X.SampleID",
  distance_method = "bray",
  view_type = c("separate", "together")
)

# Plot an exclusive unique core matrix
plot_core_matrix(
  physeq = gp, 
  group_var = "SampleType", 
  percent_samples = 0.85, 
  taxa_rank = "Phylum", 
  samplecolumn = "X.SampleID",
  log_scale = TRUE
)
```

## Getting help

If you encounter a clear bug, please file an issue with a minimal reproducible example on GitHub.
