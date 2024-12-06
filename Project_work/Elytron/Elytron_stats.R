# Function to add a file ID to each TPS entry
add_file_id <- function(tps_data, file_name) {
  modified_data <- c()
  
  # Initialize an entry buffer to store each block of data until ID line
  entry_buffer <- c()
  
  for (line in tps_data) {
    entry_buffer <- c(entry_buffer, line)
    
    # If we reach an "ID=" line, we assume the entry is complete
    if (grepl("^ID=", line)) {
      # Insert the file ID after the ID line
      entry_buffer <- c(entry_buffer, paste("Original_File_ID=", file_name, sep=""))
      
      # Append the current entry to the modified data
      modified_data <- c(modified_data, entry_buffer)
      
      # Clear the buffer for the next entry
      entry_buffer <- c()
    }
  }
  
  
  return(modified_data)
}

# List of TPS files
tps_files <- c("Aesalini.TPS", "Ceratognathini.TPS", "Chiasognathini.TPS", 
               "Lamprimini.TPS", "Lucanini.TPS", "Nicagini.TPS", "Platycerini.TPS",
               "Platyceroidini.TPS", "Streptocerini.TPS", "Syndesinae.TPS")

# Initialize an empty vector to store the modified data from all files
combined_tps_data <- c()

# Loop through each TPS file, read, modify, and append the modified data
for (tps_file in tps_files) {
  # Read the TPS file
  tps_data <- readLines(tps_file)
  
  # Extract the base file name (without ".TPS") to use as the file ID
  file_name <- gsub("\\.TPS$", "", tps_file)
  
  # Modify the TPS data by adding the file ID
  modified_data <- add_file_id(tps_data, file_name)
  
  # Append the modified data to the combined data
  combined_tps_data <- c(combined_tps_data, modified_data)
}

# Write the combined modified data to a new TPS file
writeLines(combined_tps_data, "Elytron.TPS")

cat("Modified TPS data saved to: Elytron.TPS")


##################################

# Step 1: Read the TPS file and extract landmarks
# Load necessary packages
install.packages("geomorph")
library(geomorph)
library(ggplot2)
library(dplyr)
library(viridisLite)
library(viridis)

tps_file <- "Elytron.TPS"
landmark_data <- readland.tps(tps_file, specID = "ID")  # Assuming "ID" contains the specimen IDs

# Step 2: Extract the file of origin from the TPS file
lines <- readLines(tps_file)
file_origin <- gsub("^Original_File_ID=", "", lines[grep("^Original_File_ID=", lines)])  # Extract file origin

# Ensure the length of file_origin matches the number of specimens
file_origin <- factor(file_origin)  # Convert to factor for grouping

# Create a new variable for Tribe (you may need to define this based on your data)
tribe <- file_origin  # Assuming Tribe is the same as file origin for this example

# Step 3: Perform Procrustes alignment (GPA) important file!
gpa_result <- gpagen(landmark_data)

# Step 4: Perform PCA on the aligned shapes using geomorph's gm.prcomp
pca_result <- gm.prcomp(gpa_result$coords)

# Step 5: Extract the variance explained by PC1 and PC2
explained_variance <- 100 * (pca_result$sdev^2 / sum(pca_result$sdev^2))

# Step 6: Create informative axis labels
pc1_label <- paste0("PC1 (", round(explained_variance[1], 2), "% variance)")
pc2_label <- paste0("PC2 (", round(explained_variance[2], 2), "% variance)")

# Step 7: Extract PCA scores and prepare data for plotting
pca_scores <- as.data.frame(pca_result$x)  # PCA scores from gm.prcomp
colnames(pca_scores) <- paste0("PC", 1:ncol(pca_scores))  # Ensure column names are PC1, PC2, etc.
pca_scores$Tribe <- tribe  # Add Tribe as a factor

# Step 8: Create convex hulls for each tribe
hulls <- pca_scores %>%
  group_by(Tribe) %>%
  slice(chull(PC1, PC2))  # Find the convex hull points

# Step 9: Plot the PCA results with points colored by Tribe and filled convex hulls
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Tribe)) +
  geom_point(size = 0.5) +  # Plot points with size 0.5
  geom_polygon(data = hulls, aes(fill = Tribe, group = Tribe), alpha = 0.2, palette = "Set 3") +  # Convex hulls
  labs(title = "Shape of Stag Beetle Elytra",
       x = pc1_label,  # Use the modified PC1 label
       y = pc2_label) +  # Use the modified PC2 label
  theme_minimal() +
  scale_color_brewer(palette = "Set3") +  # Choose a color-blind friendly palette
  scale_fill_brewer(palette = "Set3") +    # Matching fill colors for the convex hulls
  theme(legend.key.size = unit(3, "mm")) +  # Increase legend point size
  guides(color = guide_legend(override.aes = list(size = 3)))  # Override the legend point size

#Run GPA on my device potentially and then do relevant statistics
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Tribe)) +
  geom_point(size = 0.5) +
  geom_polygon(data = hulls, aes(fill = Tribe, group = Tribe), alpha = 0.2) +  # Convex hulls
  labs(title = "Shape of Stag Beetle Elytra",
       x = pc1_label,
       y = pc2_label) +
  theme_minimal() +
  scale_color_brewer(palette = "Set3") +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.key.size = unit(3, "mm")) +
  guides(color = guide_legend(override.aes = list(size = 3)))

# Define 'tribe' once
file_origin <- factor(file_origin)
tribe <- file_origin

# Procrustes ANOVA for overall group differences
group_test <- procD.lm(gpa_result$coords ~ tribe, iter = 999)
summary(group_test)

# Pairwise Procrustes ANOVA for specific group differences
library(RRPP)
library(geomorph)
pairwise_results <- pairwise(group_test, groups = tribe)
summary(pairwise_results)

# Adjust p-values using the Benjamini-Hochberg correction as the internet said to do
p_adjusted <- p.adjust(pairwise_results$p.value, method = "BH")
p_adjusted

# Enhanced PCA plot with confidence ellipses
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Tribe)) +
  geom_point(size = 1) +
  stat_ellipse(aes(group = Tribe), level = 0.95) +  # 95% confidence ellipses
  labs(title = "PCA of Shape by Tribe with 95% Confidence Ellipses",
       x = pc1_label,
       y = pc2_label) +
  theme_minimal() +
  scale_color_brewer(palette = "Set3") +
  theme(legend.key.size = unit(3, "mm")) +
  guides(color = guide_legend(override.aes = list(size = 3)))

# Calculate Mahalanobis distances between groups
mahalanobis_distances <- morphol.disparity(gpa_result$coords ~ tribe, groups = tribe, method = "mahalanobis", iter = 999)
summary(mahalanobis_distances)  # View Mahalanobis distances and p-values

# Calculate morphological disparity within and between tribes
disparity_results <- morphol.disparity(gpa_result$coords ~ tribe, groups = tribe, iter = 999)

# Display the summary of disparity results
summary(disparity_results)

# Run disparity analysis 
disparity_results <- morphol.disparity(gpa_result$coords ~ tribe, groups = tribe, iter = 999)

# Take disparity values for each tribe
disparity_values <- disparity_results$Procrustes.var  # Variance within each tribe
tribe_names <- levels(tribe)  # Names of the tribes

# Combine into a data frame for plotting 
disparity_df <- data.frame(Tribe = tribe_names, Disparity = disparity_values)

library(ggplot2)

# Plot disparity values
ggplot(disparity_df, aes(x = Tribe, y = Disparity, fill = Tribe)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Morphological Disparity (Shape Diversity) by Tribe",
       x = "Tribe",
       y = "Disparity (Shape Variance)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3")  # Optional color palette

pairwise_summary <- summary(pairwise_results)
print(pairwise_summary) #view pairwise summary


#HeapMap visualization 
# Step 1: Pull Pairwise results from the pairwise Procrustes ANOVAComparisons and Prepare Heatmap Data
pairwise_pvalues <- pairwise_results$p.value

# Prep tribes and heatmap data
tribes <- levels(tribe)  # Unique tribe names
heatmap_data <- expand.grid(Tribe1 = tribes, Tribe2 = tribes)  # Create all pairwise combinations
heatmap_data$PValue <- NA  # Initialize with NA values

# Fill the heatmap data with the p-values
for (i in 1:(length(tribes) - 1)) {
  for (j in (i + 1):length(tribes)) {
    pair <- paste(tribes[i], tribes[j], sep = " - ")
    if (pair %in% names(pairwise_pvalues)) {
      heatmap_data$PValue[heatmap_data$Tribe1 == tribes[i] & heatmap_data$Tribe2 == tribes[j]] <- 
        pairwise_pvalues[pair]
      heatmap_data$PValue[heatmap_data$Tribe1 == tribes[j] & heatmap_data$Tribe2 == tribes[i]] <- 
        pairwise_pvalues[pair]  # Symmetric matrix
    }
  }
}

# Need to replace NA values with 1 for non-significant comparisons!
heatmap_data$PValue[is.na(heatmap_data$PValue)] <- 1

# Add significance stars based on raw p-values into the heatmap blocks
heatmap_data$Significance <- case_when(
  heatmap_data$PValue <= 0.001 ~ "***",
  heatmap_data$PValue <= 0.01 ~ "**",
  heatmap_data$PValue <= 0.05 ~ "*",
  TRUE ~ ""  # No stars for non-significant values
)

# Step 2: Filter Significant Comparisons
significant_heatmap_data <- heatmap_data %>%
  filter(PValue <= 0.05)

# Step 3: Plot the Heatmap
library(ggplot2)

ggplot(significant_heatmap_data, aes(x = Tribe1, y = Tribe2, fill = PValue)) +
  geom_tile(color = "white") +  # Heatmap tiles with white borders
  geom_text(aes(label = Significance), color = "black", size = 3) +  # Add significance stars
  scale_fill_gradient(low = "blue", high = "red", name = "P-value", limits = c(0, 0.05)) +  # Adjust color range
  labs(title = "Significant Pairwise Procrustes ANOVA Heatmap for Elytron",
       x = "Tribe 1",
       y = "Tribe 2") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size = unit(1.5, "lines")
  )


#Clean up the heatmap

# Step 1: Extract p-values and significance stars
pairwise_pvalues <- pairwise_results$p.value

# Prepare tribes and heatmap data
tribes <- levels(tribe)  # Unique tribe names
heatmap_data <- expand.grid(Tribe1 = tribes, Tribe2 = tribes)  # Create all pairwise combinations
heatmap_data$PValue <- NA  # Initialize with NA values

# Fill the heatmap data with p-values
for (i in 1:length(tribes)) {
  for (j in 1:length(tribes)) {
    pair <- paste(tribes[i], tribes[j], sep = " - ")
    if (pair %in% names(pairwise_pvalues)) {
      heatmap_data$PValue[heatmap_data$Tribe1 == tribes[i] & heatmap_data$Tribe2 == tribes[j]] <- 
        pairwise_pvalues[pair]
      heatmap_data$PValue[heatmap_data$Tribe1 == tribes[j] & heatmap_data$Tribe2 == tribes[i]] <- 
        pairwise_pvalues[pair]  # Symmetric matrix
    }
  }
}

# Replace NA values with 1 for non-significant comparisons
heatmap_data$PValue[is.na(heatmap_data$PValue)] <- 1

# Add significance stars into heatmap blocks
heatmap_data$Significance <- case_when(
  heatmap_data$PValue <= 0.001 ~ "***",
  heatmap_data$PValue <= 0.01 ~ "**",
  heatmap_data$PValue <= 0.05 ~ "*",
  TRUE ~ ""  # No stars for non-significant values
)

# Step 2: Plot the Heatmap
ggplot(heatmap_data, aes(x = Tribe1, y = Tribe2, fill = PValue)) +
  geom_tile(color = "white") +  # Heatmap tiles with white borders
  geom_text(aes(label = Significance), color = "black", size = 3) +  # Add significance stars
  scale_fill_gradient(
    low = "white", high = "blue",  # Use blue gradient for significance
    name = "P-Value Cutoff",
    breaks = c(0.001, 0.01, 0.05, 1),
    labels = c("P ≤ 0.001 (***)", "0.001 < P ≤ 0.01 (**)", "0.01 < P ≤ 0.05 (*)", "P > 0.05 (ns)")
  ) +
  labs(
    title = "Significant Pairwise Procrustes ANOVA Heatmap for Elytron",
    x = "Tribe 1",
    y = "Tribe 2"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size = unit(1.5, "lines")
  )

