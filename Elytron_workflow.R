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

# Step 3: Perform Procrustes alignment (GPA)
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
