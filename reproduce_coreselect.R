# Reproduce Coreselect Function

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(phyloseq)
library(microViz)

# Source all R files
source("R/utils-radar.R")
source("R/utils-colors.R")
source("R/coreselect.R")

# 1. Create Mock Data (Phyloseq Object)
set.seed(123)
n_samples_per_group <- 5
groups <- c("Control", "TreatmentA", "TreatmentB")
n_samples <- n_samples_per_group * length(groups)
n_taxa <- 20

# Create Taxa Names
taxa_names <- paste0("Taxon_", LETTERS[1:n_taxa])

# OTU Table
# Initialize with zeros
otu_mat <- matrix(0, nrow = n_taxa, ncol = n_samples)
rownames(otu_mat) <- taxa_names
colnames(otu_mat) <- paste0("Sample", 1:n_samples)

# Assign groups
group_vec <- rep(groups, each = n_samples_per_group)
meta_df <- data.frame(
    SampleID = colnames(otu_mat),
    Group = group_vec,
    row.names = colnames(otu_mat)
)

# Populate OTU table with patterns
# core_taxa_control: Taxon_A, Taxon_B (Present in all Control)
# core_taxa_trtA: Taxon_B, Taxon_C (Present in all TrtA)
# core_taxa_trtB: Taxon_D (Present in all TrtB)
# Taxon_E is random noise

for (isValid in 1:n_samples) {
    grp <- group_vec[isValid]

    # Taxon A: Core in Control
    if (grp == "Control") otu_mat["Taxon_A", isValid] <- 100

    # Taxon B: Core in Control and TreatmentA
    if (grp %in% c("Control", "TreatmentA")) otu_mat["Taxon_B", isValid] <- 50

    # Taxon C: Core in TreatmentA
    if (grp == "TreatmentA") otu_mat["Taxon_C", isValid] <- 200

    # Taxon D: Core in TreatmentB
    if (grp == "TreatmentB") otu_mat["Taxon_D", isValid] <- 500

    # Taxon E: Random presence (20%)
    if (runif(1) < 0.2) otu_mat["Taxon_E", isValid] <- 10
}

physeq <- phyloseq(
    otu_table(otu_mat, taxa_are_rows = TRUE),
    sample_data(meta_df)
)

print("Phyloseq object created.")

# 2. Test coreselect
tryCatch(
    {
        print("Testing coreselect...")

        # Test with 80% prevalence (should capture the core taxa defined above)
        p <- coreselect(
            physeq = physeq,
            group_var = "Group",
            percent_samples = 80,
            abundance_threshold = 0
        )

        print("Coreselect plot generated successfully!")

        # Inspect plot data (if possible, though we can't see the object internals easily in this script without breaking encapsulation or using ggplot_build)
        # But if it ran, it's a good sign.

        if (inherits(p, "ggplot")) {
            print("Verification PASSED: Result is a ggplot object.")
        } else {
            print("Verification FAILED: Result is not a ggplot object.")
        }

        # Test with different threshold (percent > 1)
        p2 <- coreselect(
            physeq = physeq,
            group_var = "Group",
            percent_samples = 0.5, # 50%
            abundance_threshold = 10 # Higher abundance required
        )
        print("Coreselect plot (50%, >10 counts) generated successfully!")
    },
    error = function(e) {
        print(paste("Verification FAILED with error:", e$message))
        traceback()
    }
)
